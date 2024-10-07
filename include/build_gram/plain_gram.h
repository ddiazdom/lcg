//
// Created by Diaz, Diego on 3.10.2024.
//

#ifndef LCG_PLAIN_GRAM_H
#define LCG_PLAIN_GRAM_H
struct plain_gram{

    std::vector<uint64_t *> fps;
    std::vector<size_t> fps_len;
    phrase_set<uint8_t> ter_dict;
    std::vector<phrase_set<uint32_t>> nt_dicts;
    std::vector<uint32_t> comp_string;
    size_t n_levels=0;

    explicit plain_gram(size_t lvl, uint8_t sep_sym){
        assert(lvl>2);
        fps.resize(lvl);
        fps_len.resize(lvl);
        nt_dicts.resize(lvl-2);

        size_t alpha_size = std::numeric_limits<uint8_t>::max()+1;
        fps[0] = mem<uint64_t>::allocate(alpha_size);
        fps_len[0] = alpha_size;
        for(size_t i=0;i<alpha_size;i++){
            fps[0][i] = XXH3_64bits(&i, sizeof(uint8_t));
            assert(fps[0][i]!=0);
        }
        fps[0][sep_sym]=0;//small hack

        //small hack as the metasymbols are one-based
        for(size_t i=1;i<fps.size();i++){
            fps[i] = mem<uint64_t>::allocate(1);
            fps[i][0]=0;
            fps_len[i]=1;
        }
    }

    void get_gram_levels() {
        n_levels = !ter_dict.empty();
        size_t l=0;
        while(!nt_dicts[l].empty()){
            n_levels++;
            l++;
        }
    }

    void swap(plain_gram& other){
        fps.swap(other.fps);
        fps_len.swap(other.fps_len);
        ter_dict.swap(other.ter_dict);
        nt_dicts.swap(other.nt_dicts);
        comp_string.swap(other.comp_string);
        std::swap(n_levels, other.n_levels);
    }

    size_t mem_usage(){

        size_t bytes = ter_dict.mem_usage();

        for(auto const& dict: nt_dicts){
            bytes+=dict.mem_usage();
        }

        for(auto const& f_len : fps_len){
            bytes+=f_len*sizeof(uint64_t);
            bytes+=sizeof(uint64_t);//the length
        }

        return bytes;
    }

    size_t eff_mem_usage(){
        size_t bytes = ter_dict.eff_mem_usage();
        for(auto const& dict: nt_dicts){
            bytes+=dict.eff_mem_usage();
        }

        for(auto const& f_len : fps_len){
            bytes+=f_len*sizeof(uint64_t);
            bytes+=sizeof(uint64_t);//the length
        }
        return bytes;
    }

    size_t av_bytes(){
        size_t bytes = ter_dict.buff_bytes_available();
        for(auto const& dict: nt_dicts){
            bytes+=dict.buff_bytes_available();
        }
        return bytes;
    }

    inline void update_fps(size_t round){
        if(round>0){
            nt_dicts[round-1].update_fps(fps[round], fps_len[round], fps[round+1], fps_len[round+1]);
        }else{
            ter_dict.update_fps(fps[round], fps_len[round], fps[round+1], fps_len[round+1]);
        }
    }

    void clear(){
        ter_dict.clear();
        for(auto &dict: nt_dicts){
            dict.clear();
        }

        for(size_t i=1;i<fps.size();i++){
            fps[i] = mem<uint64_t>::reallocate(fps[i], 1);
            //small hack as the metasymbols are one-based
            fps[i][0]=0;
            fps_len[i]=1;
        }
    }

    void print_stats(){
        std::cout<<"Level 1, number of phrases: "<<ter_dict.size()<<",  number of symbols: "<<ter_dict.tot_symbols()<<std::endl;
        size_t round=2;
        for(auto const& lvl_set: nt_dicts){
            if(lvl_set.empty()) break;
            std::cout<<"Level "<<round<<", number of phrases: "<<lvl_set.size()<<", number of symbols: "<<lvl_set.tot_symbols()<<std::endl;
            round++;
        }
        std::cout<<"Tot. strings: "<<comp_string.size()<<std::endl;
    }

    ~plain_gram(){
        if(!ter_dict.empty()){
            std::cout<<" dict ter "<<ter_dict.size()<<", "<<ter_dict.load_factor()<<std::endl;
            ter_dict.psl_dist();
        }
        size_t i=0;
        for(auto &nt_dict : nt_dicts){
            if(!nt_dict.empty()){
                std::cout<<i<<" "<<nt_dict.size()<<", "<<nt_dict.load_factor()<<std::endl;
                nt_dict.psl_dist();
            }
            i++;
        }
    }
};
#endif //LCG_PLAIN_GRAM_H
