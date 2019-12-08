/* vcf-split.c */
int main(int argc, const char *argv[]);
void usage(const char *argv[]);
int split_vcf(const char *argv[], FILE *infile, size_t first_col, size_t last_col, size_t max_calls);
void write_output_files(const char *argv[], FILE *infile, const char *sample_ids[], size_t first_col, size_t last_col, size_t max_calls);
int split_line(const char *argv[], FILE *infile, FILE *outfiles[], const char *sample_ids[], size_t first_col, size_t last_col);
