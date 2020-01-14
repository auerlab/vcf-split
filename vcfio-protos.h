/* vcfio.c */
void vcf_skip_header(const char *argv[], FILE *infile);
void vcf_get_sample_ids(const char *argv[], FILE *infile, char *sample_ids[], size_t first_col, size_t last_col);
int vcf_read_static_fields(const char *argv[], FILE *vcf_stream, vcf_call_t *vcf_call);
int vcf_read_ss_call(const char *argv[], FILE *vcf_stream, vcf_call_t *vcf_call);
size_t vcf_read_duplicate_calls(const char *argv[], FILE *vcf_stream, vcf_duplicate_call_t *vcf_duplicate_calls);
