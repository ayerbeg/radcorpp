// plot_rc.C -- read the radcor TTree and draw the RC factor and the tail fractions of one spectrum.
// usage: root -l -b -q 'plot_rc.C("details.root", 4, "rc_root.pdf")'
void plot_rc(const char* file = "details.root", int spec = 4, const char* out = "rc_root.pdf") {
  TFile f(file);
  auto* t = f.Get<TTree>("radcor");
  TString sel = Form("spectrum==%d", spec);

  auto* c = new TCanvas("c", "radcor", 900, 380);
  c->Divide(2, 1);
  c->cd(1);
  // RC factor = measured / Born, at the kinematics actually used
  t->Draw("rc:nu", sel, "goff");
  auto* g1 = new TGraph(t->GetSelectedRows(), t->GetV2(), t->GetV1());
  g1->SetTitle(Form("spectrum %d;#nu [MeV];#sigma_{meas}/#sigma_{Born}", spec));
  g1->SetLineColor(kAzure + 2); g1->SetLineWidth(2);
  g1->Draw("AL");

  c->cd(2);
  // contributions to the radiated cross section, as fractions of the measured one
  t->Draw("tail_es/xs_in:nu", sel, "goff");
  auto* gb = new TGraph(t->GetSelectedRows(), t->GetV2(), t->GetV1());
  t->Draw("tail_ep/xs_in:nu", sel, "goff");
  auto* gc = new TGraph(t->GetSelectedRows(), t->GetV2(), t->GetV1());
  t->Draw("(rad_model-tail_es-tail_ep)/xs_in:nu", sel, "goff");
  auto* gs = new TGraph(t->GetSelectedRows(), t->GetV2(), t->GetV1());
  auto* mg = new TMultiGraph("mg", ";#nu [MeV];fraction of #sigma_{meas}");
  gs->SetLineColor(kAzure + 2); gb->SetLineColor(kOrange + 7); gc->SetLineColor(kGreen + 2);
  for (auto* g : {gs, gb, gc}) { g->SetLineWidth(2); mg->Add(g, "L"); }
  mg->Draw("A");
  auto* leg = new TLegend(0.45, 0.70, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->AddEntry(gs, "soft factor #times #tilde{#sigma} (soft+virtual)", "l");
  leg->AddEntry(gb, "tail_es (E_{s}' tail)", "l");
  leg->AddEntry(gc, "tail_ep (E_{p}' tail)", "l");
  leg->Draw();
  c->SaveAs(out);
}
