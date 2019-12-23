/* vcf-split.c */
int main(int argc, const char *argv[]);
int vcf_split(const char *argv[], FILE *vcf_infile, const char *prefix, size_t first_col, size_t last_col, size_t max_calls);
void write_output_files(const char *argv[], FILE *vcf_infile, const char *sample_ids[], const char *prefix, size_t first_col, size_t last_col, size_t max_calls);
int split_line(const char *argv[], FILE *vcf_infile, FILE *vcf_outfiles[], const char *sample_ids[], size_t first_col, size_t last_col);
id_list_t *read_selected_sample_ids(const char *argv[], const char *samples_file);
size_t read_string(FILE *fp, char *buff, size_t maxlen);
void usage(const char *argv[]);
