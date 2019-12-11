/* vcfio.c */
void vcf_skip_header(const char *argv[], FILE *infile);
void vcf_get_sample_ids(const char *argv[], FILE *infile, char *sample_ids[], size_t first_col, size_t last_col);
int vcf_read_call(const char *argv[], FILE *vcf_stream, vcf_call_t *vcf_call);
