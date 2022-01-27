#include <TFile.h>
#include <TEveManager.h>
#include <TEvePointSet.h>
#include <TEveRGBAPalette.h>
#include <TColor.h>

void pointset_retrive(TString nEvent = "0")
{
   std::cout << "Dispaying event#" << nEvent << std::endl;

   //Opening validation file
   TFile* f = new TFile("../build/src/validation.root", "read");
   TEveManager::Create();
   TString event_string("event" + nEvent);

   //Visualize clusters
   TEvePointSet* ps = new TEvePointSet();
   f->GetObject(event_string+"/Ps_clusters", ps);
   ps->SetOwnIds(kTRUE);
   gEve->AddElement(ps);
   gEve->Redraw3D();

   //Visualize hits in clusters
   TEvePointSet* ps2 = new TEvePointSet();
   f->GetObject(event_string+"/Ps_clHits", ps2);
   ps2->SetOwnIds(kTRUE);
   gEve->AddElement(ps2);
   gEve->Redraw3D();

   //Visualize outliers
   TEvePointSet* ps3 = new TEvePointSet();
   f->GetObject(event_string+"/Ps_outliers", ps3);
   ps3->SetOwnIds(kTRUE);
   gEve->AddElement(ps3);
   gEve->Redraw3D();

   f->Close();
   gEve->CloseEveWindow();

}

