import FWCore.ParameterSet.Config as cms

config = dict()

#--------- general ----------#

#--------- Set Just one to true ----------#
config["RUNONMC"] = False
#-----------------------------------------#
config["USEHAMMER"] = False
config["VERBOSE"] = False

#--------- For taus ----------#
config["DZCUT"] = 0.25 # this is fixed !!
config["FSIGCUT"] = 0
config["VPROBCUT"] = 0.005
config["TAU_CHARGE"] = 1

config["USEJSON"] = not (config["RUNONMC"])
config["JSONFILE"] = "JSON/Cert_326381-326859_HI_PromptReco_Collisions18_JSON_HF_and_MuonPhys.txt" # subdata HI 2018
#config["JSONFILE"] = "JSON/Cert_326381-327564_HI_PromptReco_Collisions18_JSON_HF_and_MuonPhys.txt" # data HI 2018
#config["JSONFILE"] = "JSON/Cert_262548-263757_PromptReco_HICollisions15_JSON.txt" # data HI 2015
config["USENOHF"] = False


#--------- basic sequences ----------#
config["DOGENPARTICLES"] = (True and config["RUNONMC"])
config["DOGENEVENT"] = (True and config["RUNONMC"])
config["DOPILEUP"] = (True and config["RUNONMC"])
config["DOVERTICES"] = True
config["DOMISSINGET"] = False

config["DOBSTAUTAU"] = True
config["DOtwoPiBSTAUTAU"] = False
config["ISTRUTH"] = False

config["DOGENHIST"] = (False and config["RUNONMC"]);

#--------- JEC ----------#

config["CORRMETONTHEFLY"] = False  # at the moment JEC available just for MC Fall17

