//
// Created by Diaz, Diego on 3.7.2023.
//

#ifndef LCG_GRAMMAR_H
#define LCG_GRAMMAR_H

struct rand_sym_t{
    size_t p;//mersene prime
    size_t a;
    size_t b;

    inline size_t operator()(size_t& sym) const {
        return (a*sym + b) % p;
    }
};

struct lc_gram_buffer_t{
    uint8_t                            sigma{}; // terminal's alphabet
    size_t                             r{}; //r: number of grammar symbols (nter + ter)
    size_t                             c{}; //c: length of the right-hand of the start symbol
    size_t                             g{}; //g: sum of the rules' right-hand sides
    size_t                             max_tsym{}; //highest terminal symbol
    uint8_t                            sep_tsym{}; //separator symbol in the collection. This is 0 if the text is a single string
    std::pair<size_t, size_t>          rl_rules={0,0}; //mark the range of run-length compressed rules

    std::string                        rules_file; // rules are concatenated in this array
    std::vector<uint8_t>               terminals;
    std::vector<rand_sym_t>            seq_breakers;

    o_file_stream<size_t>              rules_buffer;
    std::vector<long>&                 str_boundaries;//
    std::vector<size_t>                lvl_rules; // number of rules generated in every parsing round
    std::vector<size_t>                lvl_size; // number of rules generated in every parsing round



    lc_gram_buffer_t(std::string& rules_f,
                     std::vector<uint8_t>& alphabet,
                     std::vector<long>& str_bd_,
                     uint8_t sep_symbol_): sep_tsym(sep_symbol_),
                                           rules_file(rules_f),
                                           rules_buffer(rules_file, BUFFER_SIZE, std::ios::out),
                                           str_boundaries(str_bd_){

        terminals = alphabet;
        sigma = terminals.size();
        max_tsym = terminals.back();
        r = max_tsym + 1;
        g = r;

        for(size_t i=0;i<=max_tsym; i++){
            rules_buffer.push_back((i<<1UL) | 1UL);
        }
    }

    ~lc_gram_buffer_t(){
        rules_buffer.close();
    }

    //void save_to_file(std::string& output_file);
    //void load_from_file(std::string &g_file);

    /*[[nodiscard]] bool is_terminal(const size_t& id) const {
        return sym_map.find(id) != sym_map.end();
    }

    [[nodiscard]] inline bool is_rl(size_t symbol) const{
        if(rl_rules.second==0) return false;
        return symbol>=rl_rules.first && symbol < (rl_rules.first + rl_rules.second);
    }

    [[nodiscard]] inline long long int parsing_level(size_t symbol) const{
        if(symbol <= max_tsym) return 0;
        for(long long int i=0;i<int(n_p_rounds);i++){
            if(lvl_rules[i]<=symbol && symbol<lvl_rules[i+1]) return i+1;
        }
        return -1;
    }*/
};

template<class sym_type>
void insert_comp_string(lc_gram_buffer_t& gram, std::string& input_file){

    i_file_stream<sym_type> ifs(input_file, BUFFER_SIZE);

    std::vector<size_t> rules_buffer;
    basic_load_vector_from_file(gram.rules_file, rules_buffer);
    assert(rules_buffer.size()==gram.g);

    rules_buffer.reserve(rules_buffer.size()+ifs.size());
    for(size_t i=0;i<ifs.size();i++){
        size_t sym = ifs.read(i)<<1UL | (i==0);
        rules_buffer.push_back(sym);
    }
    gram.r++;
    gram.g+=ifs.size();
    gram.c=ifs.size();
    ifs.close();

    assert(gram.g==rules_buffer.size());
    basic_store_vector_to_file(gram.rules_file, rules_buffer);
}

void run_length_compress(lc_gram_buffer_t &gram_buff) {


    i_file_stream<size_t> rules(gram_buff.rules_file, BUFFER_SIZE);

    std::vector<size_t> rl_rules;
    rl_rules.reserve(gram_buff.g);

    size_t new_id = gram_buff.r;
    size_t tmp_sym;
    phrase_map_t ht;
    string_t pair(2, sym_width(rules.size()));

    size_t i=gram_buff.max_tsym+2;
    size_t run_len=1;
    bool prev_flag = rules.read(i-1) & 1UL;
    size_t prev_sym = rules.read(i-1) >> 1UL;
    size_t sym;
    bool flag;
    assert(prev_flag);

    //put the terminals
    for(size_t j=0;j<=gram_buff.max_tsym;j++){
        rl_rules.push_back(rules.read(j));
        assert(rules.read(j) & 1UL);
    }

    while(i<rules.size()) {

        sym = rules.read(i)>>1UL;
        flag =  rules.read(i) & 1UL;

        if(sym!=prev_sym || flag){

            if(run_len>1){
                pair.write(0, prev_sym);
                pair.write(1, run_len);
                auto res = ht.insert(pair.data(), pair.n_bits(), 0);

                if(res.second){
                    tmp_sym = new_id++;
                    ht.insert_value_at(res.first, tmp_sym);
                }else{
                    tmp_sym = 0;
                    ht.get_value_from(res.first, tmp_sym);
                }
            } else{
                tmp_sym = prev_sym;
            }
            //n_rules+=prev_flag;
            tmp_sym = (tmp_sym<<1UL | prev_flag);
            rl_rules.push_back(tmp_sym);

            prev_flag = flag;
            prev_sym = sym;
            run_len=0;
        }
        run_len++;
        i++;
    }

    if(run_len>1){
        pair.write(0, prev_sym);
        pair.write(1, run_len);
        auto res = ht.insert(pair.data(), pair.n_bits(), 0);
        if(res.second){
            tmp_sym = new_id;
            ht.insert_value_at(res.first, tmp_sym);
        }else{
            tmp_sym = 0;
            ht.get_value_from(res.first, tmp_sym);
        }
    } else{
        tmp_sym = prev_sym;
    }
    tmp_sym = (tmp_sym<<1UL | prev_flag);
    rl_rules.push_back(tmp_sym);

    //
    const bitstream<phrase_map_t::buff_t>& stream = ht.get_data();
    key_wrapper key_w{pair.width(), ht.description_bits(), stream};
    for(auto const& phrase : ht){
        rl_rules.push_back((key_w.read(phrase, 0)<<1UL) | 1UL);
        rl_rules.push_back(key_w.read(phrase, 1)<<1UL);
    }

    gram_buff.rl_rules.first = gram_buff.r;
    gram_buff.rl_rules.second = ht.size();
    gram_buff.r += ht.size();
    gram_buff.g  = rl_rules.size();

    std::cout<<"    Stats:"<<std::endl;
    std::cout<<"      Grammar size before:        "<<rules.size()<<std::endl;
    std::cout<<"      Grammar size after:         "<<rl_rules.size()<<std::endl;
    std::cout<<"      Number of new nonterminals: "<<ht.size()<<std::endl;
    std::cout<<"      Compression ratio:          "<<float(rl_rules.size())/float(rules.size())<<std::endl;

    basic_store_vector_to_file(gram_buff.rules_file, rl_rules);
    rules.close();
}

struct lc_gram_t {

    uint8_t                            sigma{}; // terminal's alphabet
    size_t                             r{}; //r: number of grammar symbols (nter + ter)
    size_t                             c{}; //c: length of the right-hand of the start symbol
    size_t                             g{}; //g: sum of the rules' right-hand sides
    size_t                             max_tsym{}; //highest terminal symbol
    uint8_t                            sep_tsym{}; //separator symbol in the collection. This is 0 if the text is a single string
    bool                               simplified=false;

    std::vector<uint8_t>               terminals;
    std::vector<rand_sym_t>            seq_breakers;
    std::vector<size_t>                str_boundaries;//
    std::vector<size_t>                lvl_rules;// number of rules generated in every grammar level

    int_array<size_t>                  rules;
    int_array<size_t>                  rl_ptr;
    int_array<size_t>                  rule_exp;// number of rules generated in every grammar level
    std::pair<size_t, size_t>          run_len_nt;

    explicit lc_gram_t(lc_gram_buffer_t& gram_buff){

        i_file_stream<size_t> rules_buffer(gram_buff.rules_file, BUFFER_SIZE);

        sigma = gram_buff.sigma;
        r = gram_buff.r;
        g = gram_buff.g;
        c = gram_buff.c;
        max_tsym = gram_buff.max_tsym;
        sep_tsym = gram_buff.sep_tsym;

        terminals = gram_buff.terminals;
        seq_breakers = gram_buff.seq_breakers;
        lvl_rules = gram_buff.lvl_rules;

        run_len_nt = gram_buff.rl_rules;

        rules.set_width(sym_width(r));
        rules.resize(g);

        rl_ptr.set_width(sym_width(g));
        rl_ptr.resize(r-max_tsym);

        for(size_t i=0;i<=max_tsym;i++){
            rules.write(i, rules_buffer.read(i)>>1UL);
        }

        size_t rule=0;
        for(size_t i=max_tsym+1;i<rules.size();i++){
            size_t sym = rules_buffer.read(i);
            bool first = sym & 1UL;

            rules.write(i, sym>>1UL);
            if(first){
                rl_ptr.write(rule, i);
                rule++;
            }
        }
        assert(rules.size()==g);
        assert(rule==(r-(max_tsym+1)));

        rl_ptr.write(rule, g);

        size_t offset = g-c;
        str_boundaries.resize(gram_buff.str_boundaries.size());
        size_t str=0;
        for(long & str_boundary : gram_buff.str_boundaries){
            str_boundaries[str++] = str_boundary + offset;
        }
        assert(str_boundaries[0]==offset);
        assert(str_boundaries.back()==rules.size());

        size_t acc=max_tsym+1, tmp;
        for(unsigned long & lvl_rule : lvl_rules){
            tmp = lvl_rule;
            lvl_rule = acc;
            acc+=tmp;
        }
        lvl_rules.push_back(acc);
        rules_buffer.close();
    }

    size_t serialize(std::ofstream &ofs){
        size_t written_bytes=0;
        written_bytes +=serialize_elm(ofs, sigma);
        written_bytes +=serialize_elm(ofs, r);
        written_bytes +=serialize_elm(ofs, g);
        written_bytes +=serialize_elm(ofs, c);
        written_bytes +=serialize_elm(ofs, max_tsym);
        written_bytes +=serialize_elm(ofs, sep_tsym);
        written_bytes +=serialize_elm(ofs, simplified);
        written_bytes +=serialize_elm(ofs, run_len_nt.first);
        written_bytes +=serialize_elm(ofs, run_len_nt.second);
        written_bytes +=serialize_plain_vector(ofs, terminals);
        written_bytes +=serialize_plain_vector(ofs, seq_breakers);
        written_bytes +=serialize_plain_vector(ofs, str_boundaries);
        written_bytes +=serialize_plain_vector(ofs, lvl_rules);
        written_bytes +=rules.serialize(ofs);
        written_bytes +=rl_ptr.serialize(ofs);
        written_bytes +=rule_exp.serialize(ofs);
        return written_bytes;
    }

    void load(std::ifstream &ifs){
        load_elm(ifs, sigma);
        load_elm(ifs, r);
        load_elm(ifs, g);
        load_elm(ifs, c);
        load_elm(ifs, max_tsym);
        load_elm(ifs, sep_tsym);
        load_elm(ifs, simplified);
        load_elm(ifs, run_len_nt.first);
        load_elm(ifs, run_len_nt.second);
        load_plain_vector(ifs, terminals);
        load_plain_vector(ifs, seq_breakers);
        load_plain_vector(ifs, str_boundaries);
        load_plain_vector(ifs, lvl_rules);
        rules.load(ifs);
        rl_ptr.load(ifs);
        rule_exp.load(ifs);
    }

    [[nodiscard]] inline std::pair<size_t, size_t> nt2phrase(size_t sym) const {
        assert(sym>max_tsym);
        size_t pos = sym - max_tsym-1;
        return {rl_ptr[pos], rl_ptr[pos+1]-1};
    }

    [[nodiscard]] inline bool is_terminal(size_t sym) const {
        return sym<=max_tsym;
    }

    [[nodiscard]] inline bool is_rl_sym(size_t symbol) const{
        return run_len_nt.first<=symbol && symbol<(run_len_nt.first+run_len_nt.second);
    }

    [[nodiscard]] inline size_t n_terminals() const {
        return sigma;
    }

    [[nodiscard]] inline size_t n_nonterminals() const {
        return r-max_tsym;
    }

    [[nodiscard]] inline size_t comp_str_size() const {
        return c;
    }

    [[nodiscard]] inline size_t parsing_level() const {
        return 0;
    }

    [[nodiscard]] inline size_t pos2symbol(size_t idx) const{
        assert(idx<rules.size());
        return rules.read(idx);
    }

    [[nodiscard]] inline size_t start_symbol() const {
        return r-1;
    }

    void print_stats(){
        std::cout<<"Grammar stats: "<<std::endl;
        std::cout<<"  Number of strings: "<<str_boundaries.size()-1<<std::endl;
        std::cout<<"  Number of terminals: "<<terminals.size()<<std::endl;
        std::cout<<"  Number of non-terminals: "<<n_nonterminals()<<" ("<<float(rl_ptr.size()*rl_ptr.width())/8000000<<" MBs in pointers)"<<std::endl;
        std::cout<<"  Grammar size: "<<rules.size()<<" ("<<float(rules.size()*rules.width())/8000000<<" MBs)"<<std::endl;
        std::cout<<"  Grammar breakdown: "<<std::endl;
        size_t last_nt, first_nt, tot_sym, n_rules;
        for(size_t i=0;i<lvl_rules.size()-1; i++){
            n_rules = lvl_rules[i+1]-lvl_rules[i];//number of rules in the level
            if(n_rules>0){
                first_nt = lvl_rules[i];
                last_nt = lvl_rules[i + 1] - 1;
                tot_sym = rl_ptr[last_nt - max_tsym] - rl_ptr[first_nt - (max_tsym - 1)];
                std::cout << "    Level " << (i + 1) << ": number of rules: " << n_rules << ",  number of symbols: " << tot_sym << std::endl;
            }
        }
    }

    std::pair<std::vector<uint8_t>, size_t> mark_disposable_symbols() {

        //compute which nonterminals are repeated and
        // which have a replacement of length 1
        std::vector<uint8_t> rep_nts(r + 1, 0);
        for(size_t rule=max_tsym+1;rule<r;rule++){
            auto range = nt2phrase(rule);
            if(is_rl_sym((rule))){
                rep_nts[rules[range.first]] = 2;
            }else{
                for(size_t i=range.first;i<=range.second;i++){
                    rep_nts[rules[i]]+=rep_nts[rules[i]]<2;
                }
            }
        }

        std::vector<uint8_t> rem_nts(r + 1, 0);
        //mark the rules to remove
        //1) nonterminals with freq 1
        //2) terminal symbols between [min_sym..max_sym] with
        // frequency zero: to compress the alphabet
        bool remove;
        size_t n_rem=0;
        size_t start_sym = start_symbol();
        for(size_t sym=0;sym<start_sym;sym++){
            remove = rep_nts[sym]==0 || (rep_nts[sym]==1 && sym > max_tsym && !is_rl_sym(sym));
            rem_nts[sym] = remove;
            n_rem+=remove;
        }
        return {rem_nts, n_rem};
    }

    void simplify_grammar() {

        assert(!simplified);

        auto rem_syms = mark_disposable_symbols();
        int_array<size_t> new_rules(g-rem_syms.second, sym_width(r-rem_syms.second));
        int_array<size_t> new_rl_ptrs(r-rem_syms.second, sym_width((g+1)-rem_syms.second));
        int_array<size_t> offsets(r, sym_width(r));

        size_t del_syms=0;
        for(size_t sym=0;sym<r;sym++){
            offsets[sym] = del_syms;
            del_syms+=rem_syms.first[sym];
        }

        for(size_t ter=0;ter<=max_tsym;ter++){
            if(!rem_syms.first[ter]){
                new_rules.push_back(ter-offsets[ter]);
            }
        }

        std::stack<size_t> stack;
        size_t start_sym = start_symbol();
        for(size_t sym=max_tsym+1;sym<start_sym;sym++){
            if(!rem_syms.first[sym]){

                auto range = nt2phrase(sym);
                new_rl_ptrs.push_back(new_rules.size());

                if(is_rl_sym(sym)){
                    new_rules.push_back(rules[range.first]-offsets[rules[range.first]]);
                    new_rules.push_back(rules[range.second]);
                }else{
                    for(size_t j=range.first;j<=range.second;j++){
                        if(rem_syms.first[rules[j]]) {
                            stack.push(rules[j]);
                            while(!stack.empty()){
                                auto nt = stack.top();
                                stack.pop();
                                if(rem_syms.first[nt]){
                                    assert(nt>max_tsym);
                                    assert(!is_rl_sym(nt));
                                    auto range2 = nt2phrase(nt);
                                    for(size_t k=range2.second+1;k-->range2.first;){
                                        stack.push(rules[k]);
                                    }
                                }else{
                                    new_rules.push_back(nt - offsets[nt]);
                                }
                            }
                        }else{
                            new_rules.push_back(rules[j] - offsets[rules[j]]);
                        }
                    }
                }
            }
        }

        //deal with the start symbol
        auto range = nt2phrase(start_sym);
        size_t str=0;

        new_rl_ptrs.push_back(new_rules.size());
        for(size_t j=range.first;j<=range.second;j++){
            str_boundaries[str++] = new_rules.size();
            if(rem_syms.first[rules[j]]) {
                stack.push(rules[j]);
                while(!stack.empty()){
                    auto nt = stack.top();
                    stack.pop();
                    if(rem_syms.first[nt]){
                        auto range2 = nt2phrase(nt);
                        for(size_t k=range2.second+1;k-->range2.first;){
                            stack.push(rules[k]);
                        }
                    }else{
                        new_rules.push_back(nt - offsets[nt]);
                    }
                }
            }else{
                new_rules.push_back(rules[j] - offsets[rules[j]]);
            }
        }
        new_rl_ptrs.push_back(new_rules.size());

        //TODO testing
        //size_t rm_nt=0;
        //for(size_t i=max_tsym+1;i<r;i++){
        //    if(rem_syms.first[i]) rm_nt++;
        //}
        //

        size_t del_nt = (rem_syms.second-(max_tsym+1-(sigma-1)));
        size_t del_ter = rem_syms.second - del_nt;
        std::cout<<del_ter<<std::endl;
        //std::cout<<new_rules.size()<<" "<<rules.size()<<" "<<rules.size()-new_rules.size()<<" "<<rem_syms.second<<std::endl;
        //std::cout<<new_rl_ptrs.size()<<" "<<rl_ptr.size()<<" "<<(rl_ptr.size()-new_rl_ptrs.size())<<" "<<rm_nt<<std::endl;

        assert(new_rules.size()==rules.size()-rem_syms.second);
        assert(new_rl_ptrs.size()==(rl_ptr.size()-del_nt));

        float rm_per = float(rem_syms.second)/float(r)*100;
        float comp_rat = float(new_rules.size())/float(rules.size());
        std::cout<<"    Stats:"<<std::endl;
        std::cout<<"      Grammar size before:  "<<g<<std::endl;
        std::cout<<"      Grammar size after:   "<<new_rules.size()<<std::endl;
        std::cout<<"      Deleted nonterminals: "<<rem_syms.second<<" ("<<rm_per<<"%)"<<std::endl;
        std::cout<<"      Compression ratio:    "<<comp_rat<<std::endl;

        c = new_rules.size()-str_boundaries[0];
        r -= rem_syms.second;
        g = new_rules.size();

        rules.swap(new_rules);
        rl_ptr.swap(new_rl_ptrs);

        for(auto &sym : lvl_rules){
            std::cout<<sym<<" "<<sym-offsets[sym]<<" "<<offsets[sym]<<std::endl;
            sym -= offsets[sym];
        }

        run_len_nt.first -= offsets[run_len_nt.first];

        max_tsym = terminals.size();
        simplified = true;
    }
};


void check_plain_grammar(lc_gram_t& gram, std::string& uncomp_file) {

    std::cout<<"Checking the grammar produces the exact input string"<<std::endl;
    std::cout<<"  This step is optional and for debugging purposes"<<std::endl;
    std::cout<<"  Terminals:              "<<gram.n_terminals()<<std::endl;
    std::cout<<"  Number of nonterminals: "<<gram.n_nonterminals()<<std::endl;
    std::cout<<"  Compressed string:      "<<gram.comp_str_size()<<std::endl;


    i_file_stream<uint8_t> if_stream(uncomp_file, BUFFER_SIZE);

    size_t start_symbol = gram.start_symbol();
    auto res = gram.nt2phrase(start_symbol);

    std::stack<size_t> stack;

    size_t f = res.first;
    size_t l = res.second;
    size_t idx=0;

    std::string decompression;
    size_t str=0;

    for(size_t i=f; i <= l; i++) {

        stack.emplace(gram.pos2symbol(i));
        assert(stack.size()<=if_stream.size());

        while(!stack.empty()){

            auto curr_sym = stack.top() ;
            stack.pop();

            if(gram.is_terminal(curr_sym)){
                decompression.push_back((char)curr_sym);
            }else{
                auto res2 = gram.nt2phrase(curr_sym);
                if(gram.is_rl_sym(curr_sym)){
                    assert(res2.second-res2.first+1==2);
                    size_t len = gram.pos2symbol(res2.second);
                    for(size_t j=0;j<len;j++){
                        stack.emplace(gram.pos2symbol(res2.first));
                    }
                }else{
                    for(size_t j=res2.second+1; j-->res2.first;){
                        stack.emplace(gram.pos2symbol(j));
                    }
                }
            }
        }

        for(char sym : decompression){
            if(sym!=(char)if_stream.read(idx)){
                std::cout<<(int)sym<<" "<<if_stream.read(idx)<<" "<<str<<" "<<gram.str_boundaries.size()-1<<std::endl;
            }
            assert(sym==(char)if_stream.read(idx));
            idx++;
        }
        if(gram.str_boundaries[str+1]==(i+1)){
            idx++;
            str++;
        }
        decompression.clear();
    }
    std::cout<<"\tGrammar is correct!!"<<std::endl;
}

#endif //LCG_GRAMMAR_H
