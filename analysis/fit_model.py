
__all__ = ['model', 'variables']

from ROOT import RooRealVar, RooFormulaVar, RooExponential, RooArgList, RooArgSet, RooAddPdf, RooGaussian, RooProdPdf

mass_b = RooRealVar('B_M', 'B_M', 5200, 5450)
mass_b.setRange('leftband', 5200, 5279 - 50)
mass_b.setRange('signal', 5279 - 50, 5279 + 50)
mass_b.setRange('rightband', 5279 + 50, 5450)

mass_d = RooRealVar('D~0_M', 'D~0_M', 1750, 2000)

bkg_bmass_trans = RooFormulaVar('bkg_bmass_trans', 'bkg_bmass_trans', 'B_M - 5200', RooArgList(mass_b))
bkg_dmass_trans = RooFormulaVar('bkg_dmass_trans', 'bkg_dmass_trans', 'D~0_M - 5200', RooArgList(mass_d))
prompt_dmass_trans = RooFormulaVar('prompt_bmass_trans', 'prompt_bmass_trans', 'D~0_M - 5200', RooArgList(mass_d))

# Parameters
bmass_mu = RooRealVar('sig_bmass_mu', 'sig_bmass_mu', 5280, 5200, 5400)
sig_bmass_sigma = RooRealVar('sig_bmass_sigma', 'sig_bmass_sigma', 20, 10, 100)
sig_dmass_mu = RooRealVar('sig_dmass_mu', 'sig_dmass_mu', 1870, 1750, 2000)
sig_dmass_sigma = RooRealVar('sig_dmass_sigma', 'sig_dmass_sigma', 20, 10, 100)

prompt_bmass_sigma = RooRealVar('prompt_bmass_sigma', 'prompt_bmass_sigma', 20, 10, 100)
prompt_dmass_lamb = RooRealVar('prompt_bmass_lamb', 'prompt_bmass_lamb', -1.0e-03, -1, 0)

bkg_bmass_lamb = RooRealVar('bkg_bmass_lamb', 'bkg_bmass_lamb', -1.84819e-03, -1, 0)
bkg_dmass_lamb = RooRealVar('bkg_dmass_lamb', 'bkg_dmass_lamb', -1e-3, -1, 0)

# Model
sig_bmass = RooGaussian('sig_bmass', 'sig_bmass', mass_b, bmass_mu, sig_bmass_sigma)
sig_dmass = RooGaussian('sig_dmass', 'sig_dmass', mass_d, sig_dmass_mu, sig_dmass_sigma)

bkg_bmass = RooExponential('bkg_bmass', 'bkg_bmass', bkg_bmass_trans, bkg_bmass_lamb)
bkg_dmass = RooExponential('bkg_dmass', 'bkg_dmass', bkg_dmass_trans, bkg_dmass_lamb)

prompt_bmass = RooGaussian('prompt_bmass', 'prompt_bmass', mass_b, bmass_mu, prompt_bmass_sigma)
prompt_dmass = RooExponential('prompt_dmass', 'prompt_dmass', prompt_dmass_trans, prompt_dmass_lamb)

sig = RooProdPdf('sig', 'sig', sig_bmass, sig_dmass)
bkg = RooProdPdf('bkg', 'bkg', bkg_bmass, bkg_dmass)
prompt = RooProdPdf('prompt', 'prompt', prompt_bmass, prompt_dmass)

theta1 = RooRealVar('theta1', 'theta1', 0.10, 0, 1)
theta2 = RooRealVar('theta2', 'theta2', 0.3, 0, 1)

model = RooAddPdf('model', 'model', RooArgList(sig, bkg, prompt), RooArgList(theta1, theta2))
variables = RooArgSet(mass_b, mass_d)

