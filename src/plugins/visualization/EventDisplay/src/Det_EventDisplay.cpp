//
//.cooker 3D Event Display, based on.cooker 3D Display
// modified by JCB. 

#include <Det_EventDisplay.h>

// C++ libraries
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>

#include .cookerrawtree.h"

// ROOT openGL and EVE framework headers
#include "TSystem.h"
#include "TGFrame.h"
#include "TPad.h"
#include "TGraph.h"
#include "TLine.h"
#include "TMarker.h"
#include "TGLabel.h"
#include "TG3DLine.h"
#include "TGMenu.h"
#include "TRootEmbeddedCanvas.h"
#include "TCanvas.h"
#include "TGDockableFrame.h"
#include "TGeoManager.h"
#include "TEveManager.h"
#include "TGedFrame.h"
#include "TGLObject.h"
#include "TGLViewer.h"
#include "TGLSAViewer.h"
#include "TGLEmbeddedViewer.h"
#include "TRootEmbeddedCanvas.h"
#include "TGLViewerBase.h"
#include "TGLCamera.h"
#include "TGLScenePad.h"
#include "TEveGeoNode.h"
#include "TEveViewer.h"
#include "TEveScene.h"
#include "TEvePointSet.h"
#include "TEveLine.h"
#include "TEveElement.h"
#include "TGButton.h"
#include "TExec.h"
#include "TGFrame.h"
#include "TGTextEntry.h"
#include "TGTab.h"
#include "TGIcon.h"
#include "TASImage.h"
#include "TGWindow.h"
#include "TLatex.h"
#include "TError.h"

using namespace std;

// Constructor/Destructor
Det_EventDisplay::Det_EventDisplay(TTree *in, TTree *out, TFile *inf_, TFile *outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{

  // By default, there is no geometry file input (use the default)  
  choosegeo = false;

  // By default, include the logo
  stoplogo = false;

  tablecycle=0;
  chambercycle=2;

};

Det_EventDisplay::~Det_EventDisplay()
{
};

// Instantiate the geometry
void Det_EventDisplay::GenerateDet()
{

  // Import the default geometry file or use the command line geometry
  char buffer[1024];
  if (!choosegeo) geofile = .cooker_v1"; // Default geometry file

  // Write out the file name
  snprintf(buffer,1000,"%s/.cooker/shared/gdml/%s.gdml",getenv("COOKERHOME"),geofile);

  gGeoManager=gEve->GetGeometry(buffer);
  gGeoManager->DefaultColors();
  gGeoManager->CloseGeometry();

  // Wrap the geometry for EVE, making sure to incorporate all needed levels
  etn = new TEveGeoTopNode(gGeoManager,gGeoManager->GetTopNode(),0,7);
  gEve->AddElement(etn);

  std::cout<<"\n\nConfiguring initial visualization settings.\n\nIgnore TGeoManager and TGLViewerbase info/errors...\n\n";

  // Get the pointers to all the nodes you might want to edit
  findnodes();
   
 
  // Set default colors, visibility, etc.


  // Make the world volume invisible
  gGeoManager->GetTopVolume()->SetVisibility(0);
  gGeoManager->GetTopNode()->SetVisibility(0);

  // Get ROOT to shut the front door
  gGeoManager->SetVerboseLevel(-1);

  // Draw the detector
  curcam = 0; // Start with perspective view (could change this to change default)
  gEve->FullRedraw3D(kTRUE); // Get good visual settings
  camDefaults(); // Set the default views
  gEve->FullRedraw3D(kFALSE); // Do the final redraw

};


// Command line option for choosing a non-default geometry
Long_t Det_EventDisplay::useGeo(const char* cg)
{
  choosegeo = true;
  geofile = cg;   
  return 0;
}

// Command line option for suppressing the logo to save space in the viewer
Long_t Det_EventDisplay::noLogo(bool yn)
{
  stoplogo = yn;
  return 0;
}

// Methods required for cooker

Long_t Det_EventDisplay::defineHistograms()
{
  // Define histograms (none for now, can't think of any, but if you want
  // anything just let me know)
  return 0;
};

// Main startup function
Long_t Det_EventDisplay::startup()
{
  // Suppress the ROOT warnings that are meaningless
  gErrorIgnoreLevel = 2001;

  // Default save image extension
  sprintf(ext,"png");

  // Get the run info
  runinfo = (CRTRunInfo*)getFileObject("RunInfo");

  fnamedef = "Enter filename with extension (e.g. image.pdf)";

  // Add a tab in visco
  tab=addTab("Event Display");
  
  // As long as the tab is on, go nuts with laying things out
  if (tab)
    {

      // Layout frames/viewer
      TGHorizontalFrame *rframe = new TGHorizontalFrame(tab,0,0);
      vframe = new TGVerticalFrame(rframe,0,0);
      TGVerticalFrame *logo = new TGVerticalFrame(vframe,0,0);
      //      TGHorizontalFrame *wcopt = new TGHorizontalFrame(vframe,0,0);
      TGVerticalFrame *bot = new TGVerticalFrame(tab,0,0);
      TGHorizontalFrame *camframe = new TGHorizontalFrame(bot,0,0);
      V = new TGLEmbeddedViewer(rframe);

      // Make the separator for the bottom frame
      TGHorizontal3DLine *sep = new TGHorizontal3DLine(bot);

      // White background
      V->TGLViewer::UseLightColorSet();

      // Labels
      TGLabel *viewl = new TGLabel(vframe,"Detector Display Options");
      viewl->SetTextJustify(kTextCenterX | kTextCenterY);
      viewl->SetWidth(600);

      runl = new TGLabel(vframe,"Run/Event Information");
      runl->SetTextJustify(kTextCenterX | kTextCenterY);
      runl->SetWidth(600);


      // Create the particle track color key
      bool key = false;
      if (!stoplogo)
	{
	  keyframe = new TGVerticalFrame(vframe,0,0);
	  TGHorizontalFrame * top = new TGHorizontalFrame(keyframe,0,0);
	  TGHorizontalFrame * bot = new TGHorizontalFrame(keyframe,0,0);
	  TGHorizontalFrame * bot2 = new TGHorizontalFrame(keyframe,0,0);

	  key0 = new TGLabel(keyframe,"Track Color Key");
	  key1 = new TGLabel(top,"Proton");
	  key2 = new TGLabel(top,"Positron");
	  key3 = new TGLabel(top,"Electron");
	  key4 = new TGLabel(bot,"Photon");
	  key6 = new TGLabel(bot,"Neutron");
	  key7 = new TGLabel(bot,"Mu+");
	  key8 = new TGLabel(bot2,"Mu-");
	  key9 = new TGLabel(bot2,"Pi+");
	  key10 = new TGLabel(bot2,"Pi-");
	  key5 = new TGLabel(bot2,"Other");

	  key0->SetTextJustify(kTextCenterX | kTextTop);
	  key0->SetWidth(100);
	  key1->SetTextJustify(kTextCenterX | kTextTop);
	  key1->SetWidth(100);
	  key2->SetTextJustify(kTextCenterX | kTextTop);
	  key2->SetWidth(100);
	  key3->SetTextJustify(kTextCenterX | kTextTop);
	  key3->SetWidth(100);
	  key4->SetTextJustify(kTextCenterX | kTextTop);
	  key4->SetWidth(100);
	  key5->SetTextJustify(kTextCenterX | kTextTop);
	  key5->SetWidth(100);
	  key6->SetTextJustify(kTextCenterX | kTextTop);
	  key6->SetWidth(100);
	  key7->SetTextJustify(kTextCenterX | kTextTop);
	  key7->SetWidth(100);
	  key8->SetTextJustify(kTextCenterX | kTextTop);
	  key8->SetWidth(100);
	  key9->SetTextJustify(kTextCenterX | kTextTop);
	  key9->SetWidth(100);
	  key10->SetTextJustify(kTextCenterX | kTextTop);
	  key10->SetWidth(100);


	  key0->SetTextColor((Pixel_t)0x000000);
	  key1->SetTextColor(0x0000ff);
	  key2->SetTextColor(0x999933);
	  key3->SetTextColor(0xff0000);
	  key4->SetTextColor(0x00ff00);
	  key5->SetTextColor(0xcc00ff);
	  key6->SetTextColor(0x00cdcd);
	  key7->SetTextColor(0x858585);
	  key8->SetTextColor((Pixel_t)0x000000);
	  key9->SetTextColor(0x008b8b);
	  key10->SetTextColor(0xffa500);

	  top->AddFrame(key1, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
	  top->AddFrame(key2, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
	  top->AddFrame(key3, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
	  bot->AddFrame(key4, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
	  bot->AddFrame(key6, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
	  bot->AddFrame(key7, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
	  bot2->AddFrame(key8, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
	  bot2->AddFrame(key9, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
	  bot2->AddFrame(key10, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
	  bot2->AddFrame(key5, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));

	  keyframe->AddFrame(key0, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,1));
	  keyframe->AddFrame(top, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,1));
	  keyframe->AddFrame(bot, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,1));
	  keyframe->AddFrame(bot2, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,1));

	  key = true;
	}
      char bufn[100];
      if (runinfo) sprintf(bufn,"Run Number: %d",runinfo->runNumber);
      else sprintf(bufn,"Run Number: N/A");
      runNumber = new TGLabel(vframe,bufn);
      runNumber->SetTextJustify(kTextLeft | kTextTop);
      runNumber->SetWidth(100);
      runNumber->SetMaxHeight(16);
      runNumber->ChangeOptions(runNumber->GetOptions() | kFixedHeight);
      
      // Detector Display buttons
      saveButton = new TGTextButton(camframe, "Save Image",1);
      saveButton->Connect("Clicked()","Det_EventDisplay",this,"saveImage()");
      saveName = new TGTextEntry(fnamedef,camframe);
      saveName->SetDefaultSize(370,23);
      saveName->Connect("DoubleClicked()","Det_EventDisplay",this,"Clearandblk()");
      saveName->Connect("ReturnPressed()","Det_EventDisplay",this,"saveImage()");

      // //torButton->Resize(torButton->GetSize().fWidth,(torButton->GetSize().fHeight)*0.6);
      // //torButton->Layout();

      // symbButton = new TGTextButton(vframe, "SYMBs",1);
      // symbButton->Connect("Clicked()","Det_EventDisplay",this,"symbView(=true)");

      // // gtButton = new TGTextButton(vframe, "GEM Tracker",1);
      // // gtButton->Connect("Clicked()","Det_EventDisplay",this,"gtView(=true)");

      // beamButton = new TGTextButton(vframe, "Beam Pipe",1);
      // beamButton->Connect("Clicked()","Det_EventDisplay",this,"beamView(=true)");

      // tofButton = new TGTextButton(vframe, "ToFs",1);
      // tofButton->Connect("Clicked()","Det_EventDisplay",this,"tofView(=true)");

      // targetButton = new TGTextButton(vframe, "Target Chamber and Cell",1);
      // targetButton->Connect("Clicked()","Det_EventDisplay",this,"targetView(=true)");


      // Camera buttons (consolidated to a single menu button)
      camMenu = new TGPopupMenu(gClient->GetRoot());
      camMenu->AddEntry("Perspective (Free)",1);
      camMenu->AddEntry("Orthographic Top",2);
      camMenu->AddEntry("Orthographic Rear",3);
      camMenu->AddSeparator();
      camMenu->AddEntry("Reset Main Window Cameras",4);
      camButton = new TGSplitButton(camframe,new TGHotString("Change &Camera"),camMenu,false);
      camButton->Connect("ItemClicked(Int_t)","Det_EventDisplay",this,"camSignals(Int_t)");

      // // Lumi Window Button
      // LVButton = new TGTextButton(camframe,"GEM View",1);
      // LVButton->Connect("Clicked()","Det_EventDisplay",this,"LumiView()");

      // Redraw Button
      RDButton = new TGTextButton(camframe,"Redraw Scene",1);
      RDButton->Connect("Clicked()","Det_EventDisplay",this,"RedrawReq()");

      // Layout hints  (Padding order for reference L,R,T,B)
      TGLayoutHints *stackleft = new TGLayoutHints(kLHintsCenterY | kLHintsLeft,0,2,0,0);
      TGLayoutHints *stackright = new TGLayoutHints(kLHintsCenterY | kLHintsRight,2,0,0,0);
      TGLayoutHints *stacktop = new TGLayoutHints(kLHintsCenterX | kLHintsExpandX | kLHintsTop,0,0,0,3);
      TGLayoutHints *stacktoplast = new TGLayoutHints(kLHintsCenterX | kLHintsExpandX | kLHintsTop,0,0,0,15);
      TGLayoutHints *stacktight = new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,1);
      TGLayoutHints *stacktopover = new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,1);
      TGLayoutHints *stackbottom = new TGLayoutHints(kLHintsCenterX | kLHintsBottom,0,0,2,0);
      TGLayoutHints *submenu = new TGLayoutHints(kLHintsExpandX | kLHintsCenterY | kLHintsLeft, 0,2,0,0);
      TGLayoutHints *fill = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0,0,0,0);


      
      // Instantiate the necessary EVE classes for track/hit displays
      TEveManager::Create(kFALSE);
      TEveViewer *eve_v = new TEveViewer("Eve Viewer");
      eve_v->SetGLViewer(V,V->GetFrame());
      eve_v->IncDenyDestroy();
      eve_v->AddScene(gEve->GetEventScene());
      gEve->GetViewers()->AddElement(eve_v);


      // Call the detector geometry drawer
      GenerateDet();

      
      // Layout the display

      vframe->AddFrame(logo, stacktop);
      vframe->AddFrame(runl, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 0, 5));
      vframe->AddFrame(runNumber, stacktight);

      if (key) vframe->AddFrame(keyframe, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 5, 0));
      vframe->AddFrame(viewl, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 5, 5));
      if (tableNodes.size()>0)
	{	
	  tableButton = new TGTextButton(vframe, "Table",1);
	  tableButton->Connect("Clicked()","Det_EventDisplay",this,"tableView(=true)");
	  vframe->AddFrame(tableButton,stacktop);
	}
      if (chamberNodes.size()>0)
	{	
	  chamberButton = new TGTextButton(vframe, "Chamber",1);
	  chamberButton->Connect("Clicked()","Det_EventDisplay",this,"chamberView(=true)");
	  vframe->AddFrame(chamberButton,stacktop);
	}

      vframe->MapSubwindows();
      vframe->MapWindow();

      // Make the almighty and holy.cooker logo (unless told not to at command line)
      if (!stoplogo)
	{
	  char filenameBuffer[1024];
	  strcpy(filenameBuffer,getenv("COOKERHOME"));
	  strcat(filenameBuffer,"/.cooker/shared/EventDisplay.cooker.png");
	  int width = (runNumber->GetWidth())*1.45;      
	  const TGPicture *ipic =(TGPicture *)gClient->GetPicture(filenameBuffer,width,width*480/640);
	  TGIcon .cooker = new TGIcon(logo,ipic,width,width*480/640);
	  logo->AddFrame.cooker,new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 0, 0));
	}

      // Layout the bottom bar

      camframe->AddFrame(camButton, stackleft);
      //      camframe->AddFrame(LVButton, stackleft);
      camframe->AddFrame(RDButton, stackleft);
		
      //camframe->AddFrame(LEButton, stackleft);
      camframe->AddFrame(saveButton, stackright);
      camframe->AddFrame(saveName, stackright);

      // Make the bottom bar with the separator
      bot->AddFrame(sep, new TGLayoutHints(kLHintsTop | kLHintsExpandX,0,0,3,3));
      bot->AddFrame(camframe, new TGLayoutHints(kLHintsTop | kLHintsExpandX,5,0,0,0));

      // Place the main containers
      tab->AddFrame(rframe, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,0,2,0,0));
      tab->AddFrame(bot, new TGLayoutHints(kLHintsExpandX | kLHintsBottom,0,0,0,3));
      rframe->AddFrame(vframe, new TGLayoutHints( kLHintsLeft,5,5,0,0));
      rframe->AddFrame(V->GetFrame(), new TGLayoutHints(kLHintsLeft | kLHintsExpandY|kLHintsExpandX,0,0,0,0));
      

        // Make the Event Display visco tab current
      ((TGTab *) getMemoryObject("Tab Widget"))->SetTab(2);
    
    
    }

  // Reset the ROOT warning level
  gErrorIgnoreLevel = 0;


  return 0;
};

// Menu handling

// Camera switching menu
void Det_EventDisplay::camSignals(int id)
{

  //std::cout<<"\n\nI am in the menu handler!\n\n";

  switch (id)
    {
    case 1:
      Persp();
      break;
    case 2:
      OrthoTop();
      break;
    case 3:
      OrthoRear();
      break;
    case 4:
      ResetCameras();
      break;
    }
}


// Cursor control function
void Det_EventDisplay::cursorsOff()
{
  saveName->SetState(false);
};

// Camera button functions

// Standard rotatable perspective view
void Det_EventDisplay::Persp()
{
  V->TGLViewer::SetCurrentCamera(TGLViewer::kCameraPerspXOZ);
  V->TGLViewer::SetStyle(TGLRnrCtx::kFill);
  curcam = 0;
};

// Orthographic view looking from above
void Det_EventDisplay::OrthoTop() 
{
  V->TGLViewer::SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
  V->TGLViewer::SetStyle(TGLRnrCtx::kWireFrame);
  gEve->GetEventScene()->Changed();
  gEve->FullRedraw3D(kFALSE);
  curcam = 1;
};

// Orthographic view looking down the beam line
void Det_EventDisplay::OrthoRear()
{
  V->TGLViewer::SetCurrentCamera(TGLViewer::kCameraOrthoXOZ);
  V->TGLViewer::SetStyle(TGLRnrCtx::kWireFrame);
  gEve->GetEventScene()->Changed();
  gEve->FullRedraw3D(kFALSE);
  curcam = 2;
};

void Det_EventDisplay::RedrawReq()
{
  gEve->GetEventScene()->Changed();
  gEve->FullRedraw3D(kFALSE);	
};

// Reset cameras to default views
void Det_EventDisplay::ResetCameras()
{
  // Reset the cameras and then run the initial settings
  V->TGLViewer::ResetCameras();
  camDefaults();

  // Go back to user's initial camera before reset
  switch (curcam)
    {
    case 0:
      V->TGLViewer::SetCurrentCamera(TGLViewer::kCameraPerspXOZ);
      V->TGLViewer::CurrentCamera().RotateRad(-3*M_PI/16,M_PI/4);
      break;
    case 1:
      V->TGLViewer::SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
      break;
    case 2:
      V->TGLViewer::SetCurrentCamera(TGLViewer::kCameraOrthoXOZ);
      break;   
    }
   
  // Redraw
  gEve->FullRedraw3D(kFALSE);
}

// Function to set default camera views
void Det_EventDisplay::camDefaults()
{
  // Rear view
  V->TGLViewer::SetCurrentCamera(TGLViewer::kCameraOrthoXOZ);
  V->TGLViewer::ResetCurrentCamera();
  V->TGLViewer::SetOrthoCamera(TGLViewer::kCameraOrthoXOZ,1.2,0,0,M_PI/2,M_PI);

  // Top view
  V->TGLViewer::SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
  V->TGLViewer::ResetCurrentCamera();
  V->TGLViewer::SetOrthoCamera(TGLViewer::kCameraOrthoZOY,0.55,0,0,-M_PI/2,0);

  // Perspective view
  V->TGLViewer::SetCurrentCamera(TGLViewer::kCameraPerspXOZ);
  V->TGLViewer::ResetCurrentCamera();
  V->TGLViewer::CurrentCamera().RotateRad(-3*M_PI/16,M_PI/4); // Rotate to a good view
}


// Image saving function
void Det_EventDisplay::saveImage()
{

  TString fnametemp = saveName->GetText(); // saveName is a TGTextEntry

  bool gif = fnametemp.EndsWith(".gif");
  bool gifp = fnametemp.EndsWith(".gif+");
  bool jpg = fnametemp.EndsWith(".jpg");
  bool png = fnametemp.EndsWith(".png");
  bool eps = fnametemp.EndsWith(".eps");
  bool pdf = fnametemp.EndsWith(".pdf");
  bool overwrite = fnametemp.EqualTo("File already exists! Click save again to overwrite.");

  if ((((((!gif && !gifp) && !jpg) && !png) && !eps) && !pdf) && !overwrite)
    {
      saveName->SetTextColor(0xFF0000);
      saveName->SetText("Filename must have .gif, .gif+, .jpg, .png, .eps, or .pdf extension");
    }
  else
    {

      // Load filename if not preserving for overwrite
      if (!overwrite) fname = fnametemp;

      // Check if this file name already exists
      char tdir[400];
      sprintf(tdir,"%s/%s",gSystem->WorkingDirectory(),fname.Data());
      Long_t *id=0,*size=0,*flags=0,*mt=0;
      int fexists = gSystem->GetPathInfo(tdir,id,size,flags,mt);
         
      if (fexists!=0 || overwrite) // Doesn't exist or overwrite requested
	{
	  V->TGLViewer::SavePicture(fname);
	  saveName->SetText("Image saved in current directory.");
	  saveName->SetTextColor(kBlack);

	  // Reset the default extension to the user's expressed liking
	  if (gif) sprintf(ext,"gif"); //ext = "gif";
	  else if (gifp) sprintf(ext,"gif+"); //ext ="gif+";
	  else if (jpg) sprintf(ext,"jpg"); //ext ="jpg";
	  else if (png) sprintf(ext,"png"); //ext ="png";
	  else if (eps) sprintf(ext,"eps"); //ext ="eps";
	  else if (pdf) sprintf(ext,"pdf"); //ext ="pdf";
	}
      else
	{
	  saveName->SetText("File already exists! Click save again to overwrite.");
	  saveName->SetTextColor(0xFF0000);
	}
    }

};

// Helper function for image saving text box
void Det_EventDisplay::Clearandblk()
{
  saveName->Clear();
  saveName->SetTextColor(kBlack);
  if (fname) {saveName->SetText(fname);}
};

//
// Display on/off functions for various components
//


Long_t Det_EventDisplay::finalize()
{

  return 0;
};

Long_t Det_EventDisplay::viscoEvent()
{
  return 0;
};

// Event functions

Long_t Det_EventDisplay::prepare()
{
  return 0;
};

Long_t Det_EventDisplay::execute()
{

  // cout<<"\n";

  // Clear orphaned visualization objects
  gEve->ClearOrphanage();

  // Suppress the ROOT warnings that are meaningless
  gErrorIgnoreLevel = 2001;

  // Update the image filename
  char fnbuf[100];
  sprintf(fnbuf,"run_%u_event_%d.%s",(runinfo->runNumber),(eventinfo->eventNumber)-1,ext);
  saveName->SetText(fnbuf);

  //vframe->Layout();

  gEve->GetEventScene()->Changed(); // Get all changed scene elements
  gEve->FullRedraw3D(); // Redraw the scene with changes
  return 0;
};



void printallnodes(TGeoNode *node,int level)
{
     TObjArray *k=node->GetNodes();
     std::cout<<level<<" "<<node->GetName()<<"\n";
      if (k) for (int j=0;j<k->GetEntriesFast();j++)
	       printallnodes((TGeoNode*) k->At(j),level+1);
 
}

Long_t Det_EventDisplay::findnodes()
{

  int all = 0;  // Counter to keep track of if the whole geometry is found
  //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";

  int ind = 0;
  TObjArray *l=gGeoManager->GetTopVolume()->GetNodes();
  printallnodes(gGeoManager->GetTopNode(),0);
  for (int i=0;i<l->GetEntriesFast();i++)
    {
      TString name=  l->At(i)->GetName();
      std::cout<<name<<"\n";
      
      if (name.BeginsWith("al_"))
	tableNodes.push_back((TGeoNode  *)l->At(i));
      
      if (name.BeginsWith("chamber_"))
	chamberNodes.push_back((TGeoNode  *)l->At(i));
      
}
      
    
  // Return the bit count (0 if no failure)
  return (Long_t)all;

}




// Methods required for cooker
Long_t Det_EventDisplay::cmdline(char * cmd)
{
  //add cmdline handling here

  return 0;
};

extern "C"{
  Plugin *factory(TTree *in, TTree *out, TFile *inf_, TFile *outf_, TObject *p)
  {
    return (Plugin*) new Det_EventDisplay(in,out,inf_,outf_,p);
  }
}

ClassImp(Det_EventDisplay);

// Constructor/destructor for point visualization container class

visSet::visSet()
{
  //look = new TEvePointSet();
  look = new TEveLine();
};

visSet::~visSet()
{
};
