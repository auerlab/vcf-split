/* vcf-split.c */
int main(int argc, const char *argv[]);
int vcf_split(const char *argv[], FILE *vcf_infile, const char *outfile_prefix, size_t first_col, size_t last_col, id_list_t *selected_sample_ids, size_t max_calls, flag_t flags);
void write_output_files(const char *argv[], FILE *vcf_infile, const char *sample_ids[], _Bool selected[], const char *outfile_prefix, size_t first_col, size_t last_col, size_t max_calls, flag_t flags);
int split_line(const char *argv[], FILE *vcf_infile, FILE *vcf_outfiles[], const char *sample_ids[], _Bool selected[], size_t first_col, size_t last_col, flag_t flags);
id_list_t *read_selected_sample_ids(const char *argv[], const char *samples_file);
size_t read_string(FILE *fp, char *buff, size_t maxlen);
void usage(const char *argv[]);
void tag_selected_columns(char *sample_ids[], id_list_t *selected_sample_ids, _Bool selected[], size_t first_col, size_t last_col);
int strptrcmp(const char **p1, const char **p2);
