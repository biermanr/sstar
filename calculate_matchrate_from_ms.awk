#!/usr/bin/awk -f

function calc_matchrate(){
    for(tgt_ind in tgt_inds){
        for(src_ind in src_inds){
            # Assumes ploidy = 2 in both tgt and src
            split(g[src_inds[src_ind]*2+0], src_gt_h1, "");
            split(g[src_inds[src_ind]*2+1], src_gt_h2, "");

            hap_match_pct = 0;
            for(ploidy_i = 0; ploidy_i < 2; ++ploidy_i){
                hap_site_num = 0;
                hap_shared_src_hom_site_num = 0;
                hap_shared_src_het_site_num = 0;
                split(g[tgt_inds[tgt_ind]*2+ploidy_i], tgt_gt, "");

                for(snp_i = 1; snp_i <= num_snps; ++snp_i){
                    if(tgt_gt[snp_i] == "1" || src_gt_h1[snp_i] == "1" || src_gt_h2[snp_i] == "1"){
                        ++hap_site_num;
                    }
                    if(tgt_gt[snp_i] == "1" && (src_gt_h1[snp_i] + src_gt_h2[snp_i] == 2)){
                        ++hap_shared_src_hom_site_num;
                    }
                    if(tgt_gt[snp_i] == "1" && (src_gt_h1[snp_i] + src_gt_h2[snp_i] == 1)){
                        ++hap_shared_src_het_site_num;
                    }
                }

                hap_match_src_allele_num = hap_shared_src_hom_site_num + 0.5*hap_shared_src_het_site_num;
                if(hap_site_num > 0){
                    hap_match_pct += 0.5 * (hap_match_src_allele_num / hap_site_num);
                } else {
                    hap_match_pct = "NA";
                    break;
                }
            }

            print "SIM-",sim_num,"_TGT-",tgt_inds[tgt_ind],"_SRC-",src_inds[src_ind],",",hap_match_pct;
        }
    }
}

BEGIN {
    # Unique per simulation
    num_samples = 22;
    num_snps = 25;
    tgt_ind_str = "0";
    src_ind_str = "1 2 3 4 5 6 7 8 9 10";

    # These are the number of header lines ms always adds
    num_headers = 3;
    per_sim_headers = 5;

    # One-time calculations
    headers = num_headers+per_sim_headers;
    spacing = num_samples+per_sim_headers;

    split(tgt_ind_str, tgt_inds);
    split(src_ind_str, src_inds);

    for(i in tgt_inds){ ++num_tgt_inds };
    for(i in src_inds){ ++num_src_inds };

    # Sanity checks
    # - make sure none of the tgt or src individuals are the same
    # - make sure none of the tgt or src individuals aren't indexed
    #   higher than possible given the ploidy and nsamp
    # TODO

    OFS="";
    print "sim,matchrate";
}

NR < headers {
    next;

} (NR-headers)%spacing < num_samples {
    g[(NR-headers)%spacing] = $1;

} (NR-headers)%spacing == num_samples {
    sim_num += 1;
    calc_matchrate()
    # Clear g for next simulation replicate
    delete g
} END {
    calc_matchrate()
}

