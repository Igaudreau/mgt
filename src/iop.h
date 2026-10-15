#define SZNAM 1000
#define NDAT 10

#ifdef DOUBLE
#define float double
#endif

using namespace std;

struct token {
  struct token *nxt;
  char key[SZNAM];
};

int idat, tmype;
char lprompt[SZNAM];
bool readf, echo, logon;
struct token *ftok, *oftok[NDAT];
FILE *datin[NDAT], *logout;

