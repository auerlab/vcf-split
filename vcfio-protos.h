/* vcfio.c */
void skip_header(const char *argv[], FILE *infile);
int read_field(const char *argv[], FILE *infile, char buff[], size_t buff_size);
void skip_field(const char *argv[], FILE *infile);
void get_sample_ids(const char *argv[], FILE *infile, char *sample_ids[], size_t first_col, size_t last_col);
void skip_rest_of_line(const char *argv[], FILE *infile);
