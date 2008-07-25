#ifndef ATLCONF_MISC_H
   #define ATLCONF_MISC_H

void GetDate(int *month, int *day, int *year, int *hour, int *min);
long GetInt(FILE *fpin, long Default, char *spc, char *expstr);
long GetIntRange(long Default, long Min, long Max, char *spc, char *expstr);
long GetIntVer(long Default, long Min, long Max, char *spc, char *expstr);
void GetString(FILE *fpin, char *Default, char *spc, char *expstr,
               int len, char *str0);
void GetStrVer(char *def, char *spc, char *expstr, int len, char *str);
int IsYes(char def, char *spc, char *expstr);
char GetChar(char def, char *spc, char *expstr);
int FileIsThere(char *nam);
void ATL_mprintf(int np, ...);
int GetFirstInt(char *ln);
long long GetFirstLong(char *ln);
double GetFirstDouble(char *ln);
int GetLastInt(char *ln);
long long GetLastLong(char *ln);
int fNumLines(char *fnam);
char *CmndResults(char *targ, char *cmnd);
int CmndOneLine(char *targ, char *cmnd, char *ln);
int GetIntBeforeWord(char *word, char *ln);
int GetScreenHeight();
void GetEnter(FILE *fpout);
int DisplayFile(char *fnam, FILE *fpout, int nlines);
int DisplayFile0(char *fnam, FILE *fpout);
int FoundInFile(char *fnam, char *str);
char *FindUname(char *targ);
enum ARCHFAM ProbeArchFam(char *targ);
void KillUselessSpace(char *str);
int IsBitSetInField(int *field, int bit);
void SetBitInField(int *field, int bit);
void KillUselessSpace(char *str);
char *NameWithoutPath(char *file);
void GetGccVers(char *gcc, int *comp, int *major, int *minor, int *patch);
int CompIsGcc(char *comp);
int CompIsAppleGcc(char *comp);
int CompIsMIPSpro(char *comp);
int CompIsSunWorkshop(char *comp);
int CompIsIBMXL(char *comp);
char *NewStringCopy(char *old);
char *NewAppendedString(char *old, char *app);

#define Mciswspace(C) ( (((C) > 8) && ((C) < 14)) || ((C) == 32) )
#define Mlowcase(C) ( ((C) > 64 && (C) < 91) ? (C) | 32 : (C) )

#define BADINT -777938

#endif
