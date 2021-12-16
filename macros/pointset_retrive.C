#include <TFile.h>
#include <TEveManager.h>
#include <TEvePointSet.h>
#include <TEveRGBAPalette.h>
#include <TColor.h>

void pointset_retrive(TString nEvent = "0")
{
   std::cout << "Dispaying event#" << nEvent << std::endl;

   TFile* f = new TFile("../build/src/validation_3.root", "read");
   TEveManager::Create();
   TString event_string("event" + nEvent);

   TEvePointSet* ps = new TEvePointSet();
   //f->GetObject("SummaryPlots/setpointrandom", ps);
   f->GetObject(event_string+"/Ps_clusters", ps);
   ps->SetOwnIds(kTRUE);
   gEve->AddElement(ps);
   gEve->Redraw3D();

   TEvePointSet* ps2 = new TEvePointSet();
   f->GetObject(event_string+"/Ps_clHits", ps2);
   ps2->SetOwnIds(kTRUE);
   gEve->AddElement(ps2);
   gEve->Redraw3D();

   TEvePointSet* ps3 = new TEvePointSet();
   f->GetObject(event_string+"/Ps_outliers", ps3);
   ps3->SetOwnIds(kTRUE);
   gEve->AddElement(ps3);
   gEve->Redraw3D();

   f->Close();
   gEve->CloseEveWindow();

}

