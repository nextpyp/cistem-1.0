#include <wx/wx.h>
#include <wx/app.h>
#include <wx/cmdline.h>
#include <cstdio>
#include "wx/socket.h"

#include "../../core/core_headers.h"
#include "../../core/socket_codes.h"

wxDEFINE_EVENT(wxEVT_COMMAND_MYTHREAD_LAUNCHJOB, wxThreadEvent);
wxDEFINE_EVENT(wxEVT_COMMAND_MYTHREAD_SENDINFO, wxThreadEvent);

#define SERVER_ID 100
#define GUI_SOCKET_ID 101
#define MASTER_SOCKET_ID 102

SETUP_SOCKET_CODES

class
JobControlApp : public wxAppConsole
{
	wxTimer *connection_timer;
	bool have_assigned_master;
	wxArrayString possible_gui_addresses;

	wxIPV4address active_gui_address;
	wxIPV4address my_address;
	wxIPV4address master_process_address;

	wxSocketServer *socket_server;
	wxSocketClient *gui_socket;
	wxSocketBase   *master_socket;

	int number_of_received_jobs;

	bool           gui_socket_is_busy;
	bool 		   gui_socket_is_connected;
	bool           gui_panel_is_connected;

	unsigned char   job_code[SOCKET_CODE_SIZE];

	long gui_port;
	long master_process_port;

	short int my_port;

	wxArrayString my_possible_ip_addresses;
	wxString my_port_string;
	wxString my_active_ip_address;
	wxString master_ip_address;
	wxString master_port;

	long number_of_slaves_already_connected;

	JobPackage my_job_package;

	void SetupServer();

	bool ConnectToGui();
	void OnGuiSocketEvent(wxSocketEvent& event);
	void OnMasterSocketEvent(wxSocketEvent& event);
	void OnServerEvent(wxSocketEvent& event);
	void OnConnectionTimer(wxTimerEvent& event);
	void CheckForConnections();
	void SendError(wxString error_to_send);
	void SendInfo(wxString info_to_send);

	void SendJobFinished(int job_number);
	void SendJobResult(JobResult *result_to_send);
	void SendJobResultQueue(ArrayofJobResults &queue_to_send);

	void SendAllJobsFinished(long total_timing_from_master);
	void SendNumberofConnections();

	void OnThreadLaunchJob(wxThreadEvent &event);
	void OnThreadSendInfo(wxThreadEvent &event);


	public:
		virtual bool OnInit();



};

class LaunchJobThread : public wxThread
{
	public:
    	LaunchJobThread(JobControlApp *handler, RunProfile wanted_run_profile, wxString wanted_ip_address, wxString wanted_port, const unsigned char *wanted_job_code, long wanted_actual_number_of_jobs) : wxThread(wxTHREAD_DETACHED)
		{
    		main_thread_pointer = handler;
    		current_run_profile = wanted_run_profile;
    		ip_address = wanted_ip_address;
    		port_number = wanted_port;
    		actual_number_of_jobs = wanted_actual_number_of_jobs;

    		for (int counter = 0; counter < SOCKET_CODE_SIZE; counter++)
    		{
    			job_code[counter] = wanted_job_code[counter];
    		}

		}
    	//~LaunchJobThread();
	protected:

    	JobControlApp *main_thread_pointer;
    	RunProfile current_run_profile;
    	wxString ip_address;
    	wxString port_number;
    	long actual_number_of_jobs;
    	unsigned char job_code[SOCKET_CODE_SIZE];

		void LaunchRemoteJob();
	 	void QueueInfo(wxString info_to_queue);
    	virtual ExitCode Entry();
};


wxThread::ExitCode LaunchJobThread::Entry()
{
	LaunchRemoteJob();
	return (wxThread::ExitCode)0;     // success
}


IMPLEMENT_APP(JobControlApp)


bool JobControlApp::OnInit()
{
	number_of_received_jobs = 0;

	wxPrintf("Running...\n");
	long counter;
	wxIPV4address my_address;
	wxIPV4address junk_address;
	wxString current_address;
	wxArrayString buffer_addresses;
	have_assigned_master = false;

	// set up the parameters for passing the gui address..

	static const wxCmdLineEntryDesc command_line_descriptor[] =
	{
			{ wxCMD_LINE_PARAM, "a", "address", "gui_address", wxCMD_LINE_VAL_STRING, wxCMD_LINE_OPTION_MANDATORY },
			{ wxCMD_LINE_PARAM, "p", "port", "gui_port", wxCMD_LINE_VAL_NUMBER, wxCMD_LINE_OPTION_MANDATORY },
			{ wxCMD_LINE_PARAM, "j", "job_code", "job_code", wxCMD_LINE_VAL_STRING, wxCMD_LINE_OPTION_MANDATORY },
			{ wxCMD_LINE_NONE }
	};


	wxCmdLineParser command_line_parser( command_line_descriptor, argc, argv);

	if (command_line_parser.Parse(true) != 0)
	{
		wxPrintf("\n\n");
		exit(0);
	}

	// get the address and port of the gui (should be command line options).

	wxStringTokenizer gui_ip_address_tokens(command_line_parser.GetParam(0),",");

	while(gui_ip_address_tokens.HasMoreTokens() == true)
	{
		current_address = gui_ip_address_tokens.GetNextToken();
		possible_gui_addresses.Add(current_address);
		if (junk_address.Hostname(current_address) == false)
		{
			MyDebugPrint(" Error: Address (%s) - not recognized as an IP or hostname\n\n", current_address);
			exit(-1);
		};
	}

	if (command_line_parser.GetParam(1).ToLong(&gui_port) == false)
	{
		MyDebugPrint(" Error: Port (%s) - not recognized as a port\n\n", command_line_parser.GetParam(1));
		exit(-1);
	}

	if (command_line_parser.GetParam(2).Len() != SOCKET_CODE_SIZE)
	{
		{
			MyDebugPrint(" Error: Code (%s) - is the incorrect length(%i instead of %i)\n\n", command_line_parser.GetParam(2), command_line_parser.GetParam(2).Len(), SOCKET_CODE_SIZE);
			exit(-1);
		}
	}


	// copy over job code.

	for (counter = 0; counter < SOCKET_CODE_SIZE; counter++)
	{
		job_code[counter] = command_line_parser.GetParam(2).GetChar(counter);
	}

	// initialise sockets

	wxSocketBase::Initialize();

	// Attempt to connect to the gui..

	active_gui_address.Service(gui_port);
	gui_socket = new wxSocketClient();

	// Setup the event handler and subscribe to most events

	gui_socket_is_busy = false;
	gui_socket_is_connected = false;
	gui_panel_is_connected = false;

	for (counter = 0; counter < possible_gui_addresses.GetCount(); counter++)
	{
		active_gui_address.Hostname(possible_gui_addresses.Item(counter));
		wxPrintf("\n JOB CONTROL: Trying to connect to %s:%i (timeout = 4 sec) ...\n", active_gui_address.IPAddress(), active_gui_address.Service());

		gui_socket->Connect(active_gui_address, false);
		gui_socket->WaitOnConnect(4);

		if (gui_socket->IsConnected() == false)
		{
		   gui_socket->Close();
		   wxPrintf("Connection Failed.\n\n");
		}
		else
		{
			break;
		}
	}


	if (gui_socket->IsConnected() == false)
	{
	   gui_socket->Close();
	   wxPrintf("All connections Failed! Exiting...\n");
	   return false;
	}

	gui_socket->SetFlags(wxSOCKET_WAITALL | wxSOCKET_BLOCK);
	//gui_socket->SetFlags(wxSOCKET_WAITALL);
	wxPrintf(" JOB CONTROL: Succeeded - Connection established!\n\n");
	gui_socket_is_connected = true;

	// we can use this socket to get our ip_address

	buffer_addresses = ReturnIPAddress();
	current_address = ReturnIPAddressFromSocket(gui_socket);
	if (current_address.IsEmpty() == false) my_possible_ip_addresses.Add(current_address);

	for (counter = 0; counter < buffer_addresses.GetCount(); counter++)
	{
		if (buffer_addresses.Item(counter) != current_address) my_possible_ip_addresses.Add(buffer_addresses.Item(counter));
	}

	number_of_slaves_already_connected = 0;

	// subscribe to gui events..

	Bind(wxEVT_SOCKET,wxSocketEventHandler( JobControlApp::OnGuiSocketEvent), this,  GUI_SOCKET_ID );
	gui_socket->SetEventHandler(*this, GUI_SOCKET_ID);
	gui_socket->SetNotify(wxSOCKET_CONNECTION_FLAG |wxSOCKET_INPUT_FLAG |wxSOCKET_LOST_FLAG);
	gui_socket->Notify(true);


	// Job launching event..

	Bind(wxEVT_COMMAND_MYTHREAD_LAUNCHJOB, &JobControlApp::OnThreadLaunchJob, this);
	Bind(wxEVT_COMMAND_MYTHREAD_SENDINFO, &JobControlApp::OnThreadSendInfo, this);

	// Setup the connection timer, to check for connections periodically in case the events get missed..

	Bind(wxEVT_TIMER, wxTimerEventHandler( JobControlApp::OnConnectionTimer ), this, 0);
	connection_timer = new wxTimer(this, 0);
	connection_timer->Start(5000);

	return true;

}

void LaunchJobThread::LaunchRemoteJob()
{
	long counter;
	long command_counter;
	long process_counter;
	long number_of_commands_to_run;
	long number_of_commands_run = 0;
	long number_to_run_for_this_command;

	wxIPV4address address;

	// for n processes (specified in the job package) we need to launch the specified command, along with our
	// IP address, port and job code..

	wxString executable;
	wxString execution_command;


	if(current_run_profile.controller_address == "")
	{
		executable = current_run_profile.executable_name + " " + ip_address + " " + port_number + " ";
	}
	else
	{
		executable = current_run_profile.executable_name + " " + current_run_profile.controller_address + " " + port_number + " ";
	}

	for (counter = 0; counter < SOCKET_CODE_SIZE; counter++)
	{
		executable += job_code[counter];
	}


	wxMilliSleep(2000);



	if (actual_number_of_jobs + 1 < current_run_profile.ReturnTotalJobs()) number_of_commands_to_run = actual_number_of_jobs + 1;
	else
	number_of_commands_to_run = current_run_profile.ReturnTotalJobs();

	//wxPrintf("Actual = %li, running = %li\n", actual_number_of_jobs, number_of_commands_run);


	for (command_counter = 0; command_counter <  current_run_profile.number_of_run_commands; command_counter++)
	{

		if (number_of_commands_to_run - number_of_commands_run < current_run_profile.run_commands[command_counter].number_of_copies) number_to_run_for_this_command = number_of_commands_to_run - number_of_commands_run;
		else number_to_run_for_this_command = current_run_profile.run_commands[command_counter].number_of_copies;

		execution_command =  current_run_profile.run_commands[command_counter].command_to_run;
		execution_command.Replace("$command", executable);

		execution_command += "&";

		for (process_counter = 0; process_counter < number_to_run_for_this_command; process_counter++)
		{

			wxMilliSleep( current_run_profile.run_commands[command_counter].delay_time_in_ms);

			if (process_counter == 0) QueueInfo(wxString::Format("Job Control : Executing '%s' %li times.", execution_command, number_to_run_for_this_command));

			//wxThreadEvent *test_event = new wxThreadEvent(wxEVT_COMMAND_MYTHREAD_LAUNCHJOB);
			//test_event->SetString(execution_command);

			//wxQueueEvent(main_thread_pointer, test_event);
			//wxExecute(execution_command);
			system(execution_command.ToUTF8().data());
			number_of_commands_run++;
		}

	}

	// now we wait for the connections - this is taken care of as server events..

}

void  LaunchJobThread::QueueInfo(wxString info_to_queue)
{
	wxThreadEvent *test_event = new wxThreadEvent(wxEVT_COMMAND_MYTHREAD_SENDINFO);
	test_event->SetString(info_to_queue);

	wxQueueEvent(main_thread_pointer, test_event);
}


void JobControlApp::OnThreadLaunchJob(wxThreadEvent &event)
{
	if (wxExecute(event.GetString()) == 0)
	{
		SendError(wxString::Format("Error: Failed to launch (%s)", event.GetString()));
	}
}

void JobControlApp::OnThreadSendInfo(wxThreadEvent& my_event)
{
	SendInfo(my_event.GetString());
}

void JobControlApp::SendError(wxString error_to_send)
{
//	SETUP_SOCKET_CODES

	// send the error message flag

	gui_socket->SetNotify(wxSOCKET_LOST_FLAG);

	WriteToSocket(gui_socket, socket_i_have_an_error, SOCKET_CODE_SIZE);
	SendwxStringToSocket(&error_to_send, gui_socket);

	gui_socket->SetNotify(wxSOCKET_LOST_FLAG | wxSOCKET_INPUT_FLAG);
}

void JobControlApp::SendInfo(wxString info_to_send)
{
//	SETUP_SOCKET_CODES

	// send the error message flag

	gui_socket->SetNotify(wxSOCKET_LOST_FLAG);

	WriteToSocket(gui_socket, socket_i_have_info, SOCKET_CODE_SIZE);
	SendwxStringToSocket(&info_to_send, gui_socket);

	gui_socket->SetNotify(wxSOCKET_LOST_FLAG | wxSOCKET_INPUT_FLAG);
}


void JobControlApp::SetupServer()
{
	wxIPV4address my_address;
	wxIPV4address buffer_address;

	//MyDebugPrint("Setting up Server...");

	for (short int current_port = START_PORT; current_port <= END_PORT; current_port++)
	{

		if (current_port == END_PORT)
		{
			wxPrintf("JOB CONTROL : Could not find a valid port !\n\n");
			ExitMainLoop();
			return;
		}

		my_port = current_port;
		my_address.Service(my_port);


		socket_server = new wxSocketServer(my_address, wxSOCKET_WAITALL | wxSOCKET_BLOCK);
		socket_server->SetTimeout(5);

		if (socket_server->IsOk() == true)
		{
			Bind(wxEVT_SOCKET,wxSocketEventHandler( JobControlApp::OnServerEvent), this,  SERVER_ID);

			socket_server->SetEventHandler(*this, SERVER_ID);
			socket_server->SetNotify(wxSOCKET_CONNECTION_FLAG);
			socket_server->Notify(true);

			my_port_string = wxString::Format("%hi", my_port);
			break;
		}
		else
		{
			socket_server->Destroy();
		}

	}

}

void JobControlApp::OnServerEvent(wxSocketEvent& event)
{
	  CheckForConnections();
}

void JobControlApp::CheckForConnections()
{
	wxSocketBase *sock = NULL;

	 while (1==1) // sometimes, multiple connections only seem to generate one event.. so we keep checking until there are no more connections..
	 {
		  if (socket_server == NULL) break;

		  sock = socket_server->Accept(false);

		  if (sock == NULL) break;

		  sock->SetNotify(false);
		  sock->SetFlags(wxSOCKET_WAITALL | wxSOCKET_BLOCK);
		  //sock->SetFlags(wxSOCKET_WAITALL);

		  // request identification..
		  WriteToSocket(sock, socket_please_identify, SOCKET_CODE_SIZE);
		  ReadFromSocket(sock, &socket_input_buffer, SOCKET_CODE_SIZE);

		  if ((memcmp(socket_input_buffer, job_code, SOCKET_CODE_SIZE) != 0) )
		  {

			  SendError("Unknown Job ID (Job Control), leftover from a previous job? - Closing Connection");

				  // incorrect identification - close the connection..
				  sock->Destroy();
				  sock = NULL;
		  }
		  else
		  {

			  // one of the slaves has connected to us.  If it is the first one then
			  // we need to make it the master, tell it to start a socket server
			  // and send us the address so we can pass it on to all future slaves.
			  // If we have already assigned the master, then we just need to send it
			  // the masters address.

			  if (have_assigned_master == false)  // we don't have a master, so assign it
			  {

				  master_socket = sock;
				  have_assigned_master = true;

				  WriteToSocket(sock, socket_you_are_the_master, SOCKET_CODE_SIZE);
				  master_ip_address = ReceivewxStringFromSocket(sock);
				  master_port = ReceivewxStringFromSocket(sock);

				  // setup events on the master..

				  Bind(wxEVT_SOCKET, wxSocketEventHandler( JobControlApp::OnMasterSocketEvent), this,  MASTER_SOCKET_ID);

				  sock->SetEventHandler(*this, MASTER_SOCKET_ID);
				  sock->SetNotify(wxSOCKET_INPUT_FLAG |wxSOCKET_LOST_FLAG);
				  sock->Notify(true);

				  number_of_slaves_already_connected++;
				  SendNumberofConnections();

			  }
			  else  // we have a master, tell this slave who it's master is.
			  {
				  WriteToSocket(sock, socket_you_are_a_slave, SOCKET_CODE_SIZE);
				  SendwxStringToSocket(&master_ip_address, sock);
				  SendwxStringToSocket(&master_port, sock);

				  // that should be the end of our interactions with the slave
				  // it should disconnect itself, we won't even bother
				  // setting up events for it..

				  number_of_slaves_already_connected++;
				  SendNumberofConnections();

			  }
		  }
	  }


}

void JobControlApp::OnConnectionTimer(wxTimerEvent& event)
{
	wxPrintf("Connection Timer Fired\n");
	CheckForConnections();
}

void JobControlApp::OnMasterSocketEvent(wxSocketEvent& event)
{
	//  wxString s = _("JOB CONTROL : OnSocketEvent: ");
	  wxSocketBase *sock = event.GetSocket();

	  sock->SetFlags(wxSOCKET_WAITALL | wxSOCKET_BLOCK);
	  //sock->SetFlags(wxSOCKET_WAITALL);



	  // First, print a message
//	  switch(event.GetSocketEvent())
//	  {
//	    case wxSOCKET_INPUT : s.Append(_("wxSOCKET_INPUT\n")); break;
//	    case wxSOCKET_LOST  : s.Append(_("wxSOCKET_LOST\n")); break;
//	    default             : s.Append(_("Unexpected event !\n")); break;
//	  }

	  //m_text->AppendText(s);

	  //MyDebugPrint(s);

	  // Now we process the event
	  switch(event.GetSocketEvent())
	  {
	    case wxSOCKET_INPUT:
	    {
	  	  MyDebugAssertTrue(sock == master_socket, "Master Socket event from Non Master socket??");

	      // We disable input events, so that the test doesn't trigger
	      // wxSocketEvent again.

	      sock->SetNotify(wxSOCKET_LOST_FLAG);
	      ReadFromSocket(sock, &socket_input_buffer, SOCKET_CODE_SIZE);

		  if (memcmp(socket_input_buffer, socket_send_job_details, SOCKET_CODE_SIZE) == 0)
		  {
			  // send the job details through...

			  my_job_package.SendJobPackage(sock);

		  }
		  else
		  if (memcmp(socket_input_buffer, socket_i_have_an_error, SOCKET_CODE_SIZE) == 0) // identification
		  {
			 wxString error_message;
			 error_message = ReceivewxStringFromSocket(sock);

			 // send the error message up the chain..

			 SendError(error_message);
		 }
		 else
		 if (memcmp(socket_input_buffer, socket_i_have_info, SOCKET_CODE_SIZE) == 0) // identification
		 {
			 wxString info_message;
			 info_message = ReceivewxStringFromSocket(sock);

			 // send the error message up the chain..

			 SendInfo(info_message);
		 }
		 else
		 if (memcmp(socket_input_buffer, socket_job_result, SOCKET_CODE_SIZE) == 0) // identification
		 {
			 // which job is finished and how big is the result..

			 JobResult temp_job;
			 temp_job.ReceiveFromSocket(sock);



			// if (temp_job.result_size > 0)
			 //{
				 SendJobResult(&temp_job);
			 //}

		 }
		 else
		 if (memcmp(socket_input_buffer, socket_job_result_queue, SOCKET_CODE_SIZE) == 0) // identification
		 {
			// wxPrintf("(Controller) - Received socket_job_result_queue, recieving queue\n");
			 ArrayofJobResults temp_array;

			 ReceiveResultQueueFromSocket(sock, temp_array);
			 //wxPrintf("(Controller) - Received queue of %li jobs (%i), sending on to GUI\n", temp_array.GetCount(), temp_array.Item(0).job_number);
			 SendJobResultQueue(temp_array);
			 //wxPrintf("(Controller) - Finished Sending to GUI\n");

			 number_of_received_jobs+=temp_array.GetCount();
			 //wxPrintf("Have received and sent on %i jobs\n", number_of_received_jobs);

		 }
		 else
		 if (memcmp(socket_input_buffer, socket_job_finished, SOCKET_CODE_SIZE) == 0) // identification
		 {
			 // which job is finished?

			 int finished_job;
			 ReadFromSocket(sock, &finished_job, 4);

			 // send the info to the gui

			 SendJobFinished(finished_job);
		 }
		 else
		 if (memcmp(socket_input_buffer, socket_all_jobs_finished, SOCKET_CODE_SIZE) == 0) // identification
		 {
			 long total_milliseconds_from_master;
			 ReadFromSocket(sock, &total_milliseconds_from_master, sizeof(long));
			 SendAllJobsFinished(total_milliseconds_from_master);
		 }



	      // Enable input events again.

	      sock->SetNotify(wxSOCKET_LOST_FLAG | wxSOCKET_INPUT_FLAG);
	      break;
	    }

	    case wxSOCKET_LOST:
	    {

	        if (sock == master_socket)
	        {
	        	wxPrintf("JOB CONTROL : Master Socket Disconnected!!\n");
	        	sock->Destroy();
	            ExitMainLoop();
	            abort();
	        }

	        break;
	    }

	    default:
	    {

	       	wxPrintf("weird socket communication\n");
	    	abort();

	    	break;
	    }

	  }


}

void JobControlApp::SendJobFinished(int job_number)
{
	//SETUP_SOCKET_CODES

	// get the next job..
	gui_socket->SetNotify(wxSOCKET_LOST_FLAG);

	WriteToSocket(gui_socket, socket_job_finished, SOCKET_CODE_SIZE);
	WriteToSocket(gui_socket, &job_number, 4);

	gui_socket->SetNotify(wxSOCKET_LOST_FLAG | wxSOCKET_INPUT_FLAG);
}

void JobControlApp::SendJobResult(JobResult *job_to_send)
{
	//SETUP_SOCKET_CODES

	// sendjobresultcode
	gui_socket->SetNotify(wxSOCKET_LOST_FLAG);

	WriteToSocket(gui_socket, socket_job_result, SOCKET_CODE_SIZE);
	job_to_send->SendToSocket(gui_socket);

	gui_socket->SetNotify(wxSOCKET_LOST_FLAG | wxSOCKET_INPUT_FLAG);
}

void JobControlApp::SendJobResultQueue(ArrayofJobResults &queue_to_send)
{
	gui_socket->SetNotify(wxSOCKET_LOST_FLAG);

	WriteToSocket(gui_socket, socket_job_result_queue, SOCKET_CODE_SIZE);
	SendResultQueueToSocket(gui_socket, queue_to_send);

	gui_socket->SetNotify(wxSOCKET_LOST_FLAG | wxSOCKET_INPUT_FLAG);
}


void JobControlApp::SendAllJobsFinished(long total_timing_from_master)
{
	//	SETUP_SOCKET_CODES

	// get the next job..
	gui_socket->SetNotify(wxSOCKET_LOST_FLAG);

	WriteToSocket(gui_socket, socket_all_jobs_finished, SOCKET_CODE_SIZE);
	WriteToSocket(gui_socket, &total_timing_from_master, sizeof(long));

	gui_socket->SetNotify(wxSOCKET_LOST_FLAG | wxSOCKET_INPUT_FLAG);

}

void JobControlApp::SendNumberofConnections()
{
	//SETUP_SOCKET_CODES

	// get the next job..
	gui_socket->SetNotify(wxSOCKET_LOST_FLAG);

	WriteToSocket(gui_socket, socket_number_of_connections, SOCKET_CODE_SIZE);
	WriteToSocket(gui_socket, &number_of_slaves_already_connected, 4);

	gui_socket->SetNotify(wxSOCKET_LOST_FLAG | wxSOCKET_INPUT_FLAG);

	int number_of_commands_to_run;

	if (my_job_package.number_of_jobs + 1 < my_job_package.my_profile.ReturnTotalJobs()) number_of_commands_to_run = my_job_package.number_of_jobs + 1;
	else number_of_commands_to_run = my_job_package.my_profile.ReturnTotalJobs();

	if (number_of_slaves_already_connected == number_of_commands_to_run)
	{

		//wxPrintf("All connections completed\n");
		connection_timer->Stop();
		Unbind(wxEVT_TIMER, wxTimerEventHandler( JobControlApp::OnConnectionTimer ), this);
		delete connection_timer;
		Unbind(wxEVT_SOCKET,wxSocketEventHandler( JobControlApp::OnServerEvent), this,  SERVER_ID);
		socket_server->Destroy();
		socket_server=NULL;
		wxPrintf("Socket Server is now NULL\n");

	}
}

void JobControlApp::OnGuiSocketEvent(wxSocketEvent& event)
{
	  wxString s = _("JOB CONTROL : OnSocketEvent: ");
	  wxSocketBase *sock = event.GetSocket();

	  sock->SetFlags(wxSOCKET_WAITALL | wxSOCKET_BLOCK);
	  //sock->SetFlags(wxSOCKET_WAITALL);

	  MyDebugAssertTrue(sock == gui_socket, "GUI Socket event from Non GUI socket??");

	  // First, print a message
	  switch(event.GetSocketEvent())
	  {
	    case wxSOCKET_INPUT : s.Append(_("wxSOCKET_INPUT\n")); break;
	    case wxSOCKET_LOST  : s.Append(_("wxSOCKET_LOST\n")); break;
	    default             : s.Append(_("Unexpected event !\n")); break;
	  }

	  //m_text->AppendText(s);

	  //MyDebugPrint(s);

	  // Now we process the event
	  switch(event.GetSocketEvent())
	  {
	    case wxSOCKET_INPUT:
	    {
	      // We disable input events, so that the test doesn't trigger
	      // wxSocketEvent again.
	      sock->SetNotify(wxSOCKET_LOST_FLAG);
	      ReadFromSocket(sock, socket_input_buffer, SOCKET_CODE_SIZE);


	      if (memcmp(socket_input_buffer, socket_please_identify, SOCKET_CODE_SIZE) == 0) // identification
	      {
	    	  // send the job identification to complete the connection
	    	  WriteToSocket(sock, job_code, SOCKET_CODE_SIZE);
	      }
	      else
	      if (memcmp(socket_input_buffer, socket_you_are_connected, SOCKET_CODE_SIZE) == 0) // we are connected to the relevant gui panel.
	      {
	    	  gui_panel_is_connected = true;

	    	  // ask the panel to send job details..

	    	  WriteToSocket(sock, socket_send_job_details, SOCKET_CODE_SIZE);


	      }
	      else
		  if (memcmp(socket_input_buffer, socket_ready_to_send_job_package, SOCKET_CODE_SIZE) == 0) // we are connected to the relevant gui panel.
		  {
			  // receive the job details..

			  my_job_package.ReceiveJobPackage(sock);

			  wxString ip_address_string;

			  for (int counter = 0; counter < my_possible_ip_addresses.GetCount(); counter++)
			  {
			  		if (counter != 0) ip_address_string += ",";
			  		ip_address_string += my_possible_ip_addresses.Item(counter);
			  }

			  // start the server..
			  SetupServer();

			  LaunchJobThread *launch_thread = new LaunchJobThread(this, my_job_package.my_profile, ip_address_string, my_port_string, job_code, my_job_package.number_of_jobs);

			  if ( launch_thread->Run() != wxTHREAD_NO_ERROR )
			  {
				  MyPrintWithDetails("Can't create the launch thread!");
				  delete launch_thread;
				  ExitMainLoop();
				  return;
			  }
		  }
	      else
		  if (memcmp(socket_input_buffer, socket_time_to_die, SOCKET_CODE_SIZE) == 0)
		  {
			  // pass message on to master if we have one..

			  if (have_assigned_master == true)
			  {
				  master_socket->Notify(false);
				  WriteToSocket(master_socket, socket_time_to_die, SOCKET_CODE_SIZE);
				  master_socket->Destroy();
			  }

			  // destroy the server..

			  if (socket_server != NULL) socket_server->Destroy();

			  // close Gui connection..

			  sock->Destroy();

			  // exit..

			  ExitMainLoop();
			  exit(0);
			  return;
		  }


	      // Enable input events again.

	      sock->SetNotify(wxSOCKET_LOST_FLAG | wxSOCKET_INPUT_FLAG);
	      break;
	    }

	    case wxSOCKET_LOST:
	    {

	    	if (have_assigned_master == true)
	    	{
	    		master_socket->Notify(false);
	    		WriteToSocket(master_socket, socket_time_to_die, SOCKET_CODE_SIZE);
	    		 master_socket->Destroy();
	    	}

	    	// destroy the server..

	    	if (socket_server != NULL) socket_server->Destroy();

	        wxPrintf("JOB CONTROL : GUI Socket Disconnected!!\n");
	        sock->Destroy();
	        ExitMainLoop();
	        exit(0);

	        break;
	    }
	    default: ;
	  }
}


