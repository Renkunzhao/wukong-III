#include "TGraph.h"
#include "TMultiGraph.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"

int main(int argc, char **argv)
{
  TApplication app("app", &argc, argv);

  TCanvas* c = new TCanvas("c", "Something", 0, 0, 800, 600);
  auto gr = new TGraph();
  for (int i=0; i<20; i++) gr->AddPoint(i*0.1, 10*sin(i*0.1+0.2));
  gr->Draw("AL*");   // Draw() specifies the drawing option.

  TCanvas* c1 = new TCanvas("c1", "Something", 0, 0, 800, 600);
  TMultiGraph *mg = new TMultiGraph();
  auto gr1 = new TGraph();
  for (int i=0; i<20; i++) gr1->AddPoint(i*0.1, 10*cos(i*0.1+0.2));
  auto gr2 = new TGraph();
  for (int i=0; i<20; i++) gr2->AddPoint(i*0.1, 10*sin(i*0.1+0.2));

  for (int i=20; i<40; i++) gr->AddPoint(i*0.1, 10*sin(i*0.1+0.2));
  gr->Draw("AL*");   // Draw() specifies the drawing option.
  mg->Add(gr1);
  mg->Add(gr2);
  mg->Draw("AL*");

  // c->Modified(); c->Update();
  app.Run();
  return 0;
}