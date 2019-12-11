/* tsvio.c */
int tsv_read_field(const char *argv[], FILE *infile, char buff[], size_t buff_size);
void tsv_skip_field(const char *argv[], FILE *infile);
void tsv_skip_rest_of_line(const char *argv[], FILE *infile);
