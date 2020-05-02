#include "Config.h"

ClassImp(Config)

Config::Config() {
  std::ifstream configFile("config.cf");
  if(configFile.is_open()) {
    std::string line;
    while(getline(configFile,line)){

      // remove whitespace
      line.erase( std::remove_if(line.begin(), line.end(), isspace),line.end());
      
      // ignore comments
      if(line[0]=='#' || line.empty()) continue;

      auto delimPos = line.find("=");
      auto varName = line.substr(0,delimPos);
      auto varVal = line.substr(delimPos+1);

      if(varName=="Experiment") Experiment = (TString) varVal;
      if(varName=="EbeamEn") EbeamEn = std::stod(varVal);
      if(varName=="PbeamEn") PbeamEn = std::stod(varVal);
    }

  }
  else {
    printf("Config: using default CLAS settings\n");
    Experiment = "clas";
    EbeamEn = 10.6041;
    PbeamEn = 0;
  }

  // set kinematic ranges (array index: 0=min, 1=max)
  if(this->Experiment=="clas") {
    bdQ2[0] = 0;          bdQ2[1] = 12;
    bdW[0] = 0;           bdW[1] = 6;
    bdMh[0] = 0;          bdMh[1] = 3;
    bdMmiss[0] = 0;       bdMmiss[1] = 6;
    bdEta[0] = -1;        bdEta[1] = 5;
    bdHadP[0] = 0;        bdHadP[1] = 10;
  }
  else if(this->Experiment=="eic") {
    bdQ2[0] = 0;          bdQ2[1] = 200;
    bdW[0] = 0;           bdW[1] = 60;
    bdMh[0] = 0;          bdMh[1] = 5;
    bdMmiss[0] = 0;       bdMmiss[1] = 60;
    bdEta[0] = -5;        bdEta[1] = 5;
    bdHadP[0] = 0;        bdHadP[1] = 70;
  };

}


void Config::PrintVars() {
  printf("Experiment = %s\n",Experiment.Data());
  printf("EbeamEn = %f\n",EbeamEn);
  printf("PbeamEn = %f\n",PbeamEn);
};
