/* Several functions*/

void stopcode (void);
void check_stop (void);
void check_params (int argc, char *argv[]);
void splitdomain (void);
void initvars (void);
void finalize (void);
void selectposition (void);
void setbordersize (void);
void selectborder (int direction);	/*0 x-, 1 x+, 2 y-, 3 y+, 4 z-, 5 z+ */
int checkgrid (void);
void computegrid (void);
