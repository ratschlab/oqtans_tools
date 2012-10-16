char BUILD_REVERSE_INDEX;
char VERBOSE;
char HAS_SLOT;
unsigned int SLOT;
unsigned int SLOT_REV;
int INDEX_DEPTH;
long int POWER[MAX_INDEX_DEPTH];
long int BINARY_CODE[5];
int debug;
char GENOME_FILE_NAME[500];
char CHR_INDEX_FILE_NAME[500];
char MAPFWD_INDEX_FILE_NAME[500];
char MAPREV_INDEX_FILE_NAME[500];
char META_INDEX_FILE_NAME[500];
char GENOME_OUT_FILE_NAME[500];
FILE *CHR_INDEX_FP;
FILE *MAPFWD_INDEX_FP;
FILE *MAPREV_INDEX_FP;
FILE *META_INDEX_FP;
FILE *GENOME_FP;
FILE *GENOME_OUT_FP;
unsigned int CHR_LENGTH;
char *CHR_SEQ;
char CHR_DESC[2000];
char CHR_DESC_TMP[2000];
unsigned int SLOT_COUNTER; //counts the number of slots used by the forward index (equals the number of reverse slots though the slots themselves are different)
unsigned long int POSITION_COUNTER;
POS *BLOCK_TABLE;
unsigned int BLOCK;
unsigned short int POSITION;
BIN **INDEX;//[INDEX_SIZE];
//BIN **INDEX_REV;//[INDEX_SIZE];
long int NUM_USED_SLOTS; //different to SLOT_COUNTER! This counts the number of different used slots, used by reverse and forward Index.
int *USED_SLOTS;//[INDEX_SIZE];
STORAGE *MEM_MGR;

unsigned int LONGEST_CHROMOSOME;

