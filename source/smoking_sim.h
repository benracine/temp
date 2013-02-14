// CISNET (www.cisnet.cancer.gov)
// Lung Cancer Base Case Group
// Smoking History Simulation Application
// Application to Simulate Initiation and Cessation Ages of individuals based on sex, race and year of birth.
// File: smoking_sim.h
// Author: Martin Krapcho & Ben Racine
// E-Mail: KrapchoM@imsweb.com & ben.racine@cornerstonenw.com
// NCI Contacts: Rocky Feuer
// Version 6.2.3

#ifndef _SMOKING_SIM_H
#define _SMOKING_SIM_H

#include "mersenne_class.h"
#include "sim_exception.h"
#include <string.h>
#include <iostream>

// Constants used in Excess Risk Former Smokers' formula
#define B0 -0.1711
#define B1 0.00102
#define B2 0.00171
#define B3 1.08

extern short wSIM_CUTOFF_YEAR;
extern const short wMIN_IMMEDIATE_CESSATION_YEAR;
extern const char sSEX_LABELS[2][7];
extern const char sRACE_LABELS[2][10];


class Smoking_Simulator {

   // Labels and Enumerated Data Types for the class
   public:

      enum DataType {DATA_Initiation = 1, DATA_Cessation};
      enum OutputType {OUT_DataOnly = 1, OUT_TextReport, OUT_TimeLine, OUT_XML_Tags, OUT_Uninitialized};

      // Individuals smoking status
      enum SmokingStatus {SMKST_Never = 0, SMKST_Current, SMKST_Former, SMKST_NumValues};

      // Individuals smoking frequency quintile (light to heavy)
      enum SmokingIntensity {SMKR_Light = 0, SMKR_LgtMed, SMKR_Medium, SMKR_MedHvy, SMKR_Heavy, SMKR_NumGroups, SMKR_Uninitialized};

      // Columns of data in the other COD Life Table file
      enum LifeTableColumns {COL_Never = 0, COL_Current_Q1, COL_Current_Q2, COL_Current_Q3, COL_Current_Q4, COL_Current_Q5, COL_NumColumns};

      // These 2 enums are used to write the input tag for the web version output
      enum Sex {SEX_Male = 0, SEX_Female, NUM_SEXES};
      enum Race {RACE_AllRaces = 0, NUM_RACES};

 	// Private Member Variables
   private:

      // Psuedo Random Number Generator Variables
      MersenneTwister *gpInitiationPRNG;  // PRNG for Initiation Probabilities
      MersenneTwister *gpCessationPRNG;   // PRNG for Cessation Probabilities
      MersenneTwister *gpLifeTablePRNG;   // PRNG for Other COD Probabilities
      MersenneTwister *gpIndivRndsPRNG;   // PRNG for all "random" variables that only
                                          // need one value per individual (ie smoking intensity quintile)
                                          // Program will allow 20 of these values per individual

      // Probability Arrays
      double *gdInitiationProbs;  // Prob of initiation by race/sex/year of birth and age
      double *gdCessationProbs;   // Prob of cessation by race/sex/year of birth and age
      double *gdLifeTableProbs;   // Prob of COD other than lung cancer by race/sex/year of birth/age and smoking status
      double *gdIntensityProbs;   // Prob of being a light to heavy smoker (for individuals that begin smoking)

      // Cigarettes per day by race, sex, YOB and age (and smoking intensity? %bjr)
      long double *gdCigarettesPerDay;
      

      // Data limit variables
      short gwNumBirthCohorts;    // Number of birth cohorts Available
      short *gwYOBCohortStartYrs; // Starting year for each of the birth cohort groups
      short *gwYOBCohortEndYrs;   // Ending year for each of the birth cohort groups
      short gwNumRaceValues;      // Number of Races Available
      short gwNumSexValues;       // Number of Sexes Available
      short gwMinInitiationAge;   // Min initiation age (assumed constant for all cohort groups)
      short gwMinCessationAge;    // Min cessation age (assumed constant for all cohort groups)
      short gwMaxInitiationAge;   // Max initiation age (Max possible age, data may quit before max age)
      short gwMaxCessationAge;    // Max cessation age (Max possible age, data may quit before max age)
      short gwMinLifeTableAge;
      short gwMaxLifeTableAge;
      short gwMinLifeTableYear;
      short gwMaxLifeTableYear;
      short gwNumIntensityGrps;   // Number of CPD Intesity groups
      short gwIntensityMinAge;    // Minimum age among the smoking intensity group probabilities
      short gwIntensityMaxAge;    // Maximum age among the smoking intensity group probabilities
      short gwCpdMinAge;          // Minimum age in the cigarettes per day data
      short gwCpdMaxAge;          // Maximum age in the cigarettes per day data
      short gwImmediateCessYear;  // Year when all smokers automatically quit smoking. 0 = option not used.
      bool gbImmediateCessation;  // Is immediatte Cessation turned on

      // Person Variables, Store the results for the last person simulated
      short gwPersonsYOB;          // Year Of Birth
      short gwPersonsRace;         // Race
      short gwPersonsSex;          // Sex
      short gwPersonsInitAge;      // Age of Smoking Initiation
      short gwPersonsCessAge;      // Age of Smoking Cessation
      short gwPersonsAgeAtDeath;   // Age at death from COD other than lung cancer
      SmokingIntensity     gwPersonsSmkIntensity; // The smoking intesity group for the person (smokers only)
      double *gdPersonsCPDbyAge;   // Cigarettes smoked per day by age
      double gdPersonsAvgCPD;      // Average num of Cigarettes smoked per day (used for COD in former smokers)

      // Offset values for Probability Arrays
      short gwInitProbRaceOffset; // Initiation Array - Race Offset
      short gwInitProbSexOffset;  // Initiation Array - Sex Offset
      short gwInitProbYOBOffset;  // Initiation Array - YOB Offset
      short gwCessProbRaceOffset; // Cessation Array  - Race Offset
      short gwCessProbSexOffset;  // Cessation Array  - Sex Offset
      short gwCessProbYOBOffset;  // Cessation Array  - YOB Offset
      long  glLifeTabAgeOffset;   // Other COD Array  - Age Offset
      long  glLifeTabRaceOffset;  // Other COD Array  - Race Offset
      long  glLifeTabSexOffset;   // Other COD Array  - Sex Offset
      long  glLifeTabYOBOffset;   // Other COD Array  - YOB Offset
      long  gwIntensityAgeOffset; // Smoking Intesity Array - Age Offset
      long  gwIntensitySexOffset;
      long  gwIntensityRaceOffset;
      long  glCpdAgeOffset;       // Cigarettes per Day Array  - Age Offset
      long  glCpdRaceOffset;      // Cigarettes per Day Array  - Race Offset
      long  glCpdSexOffset;       // Cigarettes per Day Array  - Sex Offset
      long  glCpdYOBOffset;       // Cigarettes per Day Array  - YOB Offset

      short gwNumSmokingGrps;

      OutputType           geOutputType;

      double      gdTempIntensityProb; // Persons intensity prob, remove from final

      void Init();
      void Free();
      void CalcCigarettesPerDay();
      void CalcCigarettesPerDaySwitch();
      short GetAgeOfDeathFromOtherCOD(short wStartAge, short wEndAge, SmokingStatus eStatus, bool &bWentPastData);
      double GetNextInitRand();
      double GetNextCessRand();
      double GetNextLifeTabRand();
      double GetNextRandForIndiv();
      void InitPRNGs(unsigned long ulInitSeed, unsigned long ulCessSeed, unsigned long ulLifeTabSeed, unsigned long ulIndRndsSeed);
      void LoadCPDIntensityProbs(const char* sDataFileName);
      void LoadCPDFile(const char* sCpdDataFile);
      void LoadOtherCODFile(const char* sLifeTableFileName);
      void LoadProbabilityData(const char* sDataFileName, DataType eFileType);
      void OversamplePRNGs();

   public:
      Smoking_Simulator(const char* sInitiationProbFile, const char* sCessationProbFile,
                        const char* sLifeTableFile,      const char* sCpdIntensityProbFile,
                        const char* sCpdDataFile,        unsigned long ulInitPRNGSeed,
                        unsigned long ulCessPRNGSeed,    unsigned long ulLifeTabSeed,
                        unsigned long ulIndivRndsSeed,   short wOutputType,
                        short wCessationYear);

      ~Smoking_Simulator();

      short GetMaxYearOfBirth();
      short GetMinYearOfBirth();
      short GetNumRaceValues() { return gwNumRaceValues;};
      short GetNumSexValues() { return gwNumSexValues;};
      short GetYOBCohortGroup(short wYearBirth);

      void RunSimulation(const char* sInputFileName, const char* sOutputFileName = 0, bool bPrintToScreen = true);
      void RunSimulation(short wRace, short wSex, short wYearBirth, FILE* pOutStream = 0);

      void SetOutputType(short wOutputType);
      void WriteAsData(FILE *pOutStream);
      void WriteAsText(FILE *pOutStream);
      void WriteAsTimeline(FILE *pOutStream);
      void WriteAsXML(FILE *pOutStream);
      void WriteToStream(FILE *pOutStream);
};

#endif

