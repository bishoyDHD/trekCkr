#include "Det_EventDisplay.h"
#include "TGeoNode.h"
#include "TEveManager.h"
#include "TEveScene.h"

void Det_EventDisplay::tableView(bool rdraw)
{
  for (auto it:tableNodes)
    {
      it->SetVisibility(tablecycle);
      it->VisibleDaughters(tablecycle);
    }
  tablecycle=(tablecycle+1)%2;
  
  if (rdraw)
    {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
    }
}


void Det_EventDisplay::chamberView(bool rdraw)
{

  switch (chambercycle)
    {
    case 0:
    case 1:
      for (auto it:chamberNodes)
	{
	  it->SetVisibility(chambercycle);
	  it->VisibleDaughters(chambercycle);
	  it->GetVolume()->SetTransparency(0);
	}
    break;
    case 2:
       for (auto it:chamberNodes)
	{
	  it->SetVisibility(1);
	  it->VisibleDaughters(1);
	  it->GetVolume()->SetTransparency(75);
	}
      break;
    }
      
  chambercycle=(chambercycle+1)%3;
  
  if (rdraw)
    {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
    }
}
