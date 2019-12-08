/* vcf-split.c */
int main(int argc, const char *argv[]);
void usage(const char *argv[]);
int split_vcf(const char *argv[], FILE *infile, size_t first_col, size_t last_col, size_t max_calls);
void skip_header(const char *argv[], FILE *infile);
void get_sample_ids(const char *argv[], FILE *infile, char *sample_ids[], size_t first_col, size_t last_col);
void write_output_files(const char *argv[], FILE *infile, const char *sample_ids[], size_t first_col, size_t last_col, size_t max_calls);
void skip_rest_of_line(const char *argv[], FILE *infile);
int split_line(const char *argv[], FILE *infile, FILE *outfiles[], const char *sample_ids[], size_t first_col, size_t last_col);
int read_field(const char *argv[], FILE *infile, char buff[], size_t buff_size);
void skip_field(const char *argv[], FILE *infile);
