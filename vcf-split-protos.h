/* vcf-split.c */
int main(int argc, char *argv[]);
int vcf_split(char *argv[], FILE *vcf_infile, const char *outfile_prefix, size_t first_col, size_t last_col, id_list_t *selected_sample_ids, size_t max_calls, flag_t flags, vcf_field_mask_t field_mask);
void write_output_files(char *argv[], FILE *vcf_infile, FILE *header, const char *all_sample_ids[], _Bool selected[], const char *outfile_prefix, size_t first_col, size_t last_col, size_t max_calls, flag_t flags, vcf_field_mask_t field_mask);
int xt_split_line(char *argv[], FILE *vcf_infile, FILE *vcf_outfiles[], const char *all_sample_ids[], _Bool selected[], size_t first_col, size_t last_col, size_t max_calls, flag_t flags, vcf_field_mask_t field_mask);
void dump_line(char *argv[], const char *message, bl_vcf_t *vcf_call, size_t line_count, size_t col, size_t first_col, const char *all_sample_ids[], char *genotype);
id_list_t *read_selected_sample_ids(char *argv[], const char *samples_file);
size_t read_string(FILE *fp, char *buff, size_t maxlen);
void usage(char *argv[]);
void tag_selected_columns(char *all_sample_ids[], id_list_t *selected_sample_ids, _Bool selected[], size_t first_col, size_t last_col);
