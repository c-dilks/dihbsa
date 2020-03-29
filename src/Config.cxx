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
}


void Config::PrintVars() {
  printf("Experiment = %s\n",Experiment.Data());
  printf("EbeamEn = %f\n",EbeamEn);
  printf("PbeamEn = %f\n",PbeamEn);
};
