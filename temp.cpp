//
//  new_main.cpp
//  Polymer_cpp
//
//  Created by Baljyot Parmar on 2018-05-15.
//  Copyright © 2018 Baljyot Parmar. All rights reserved.
//
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <fcntl.h>
#include <sys/mman.h>
#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/multi_array.hpp>
#include <vector>
#include <boost/range/irange.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>
#include<boost/range/numeric.hpp>
#include <unordered_map>
#include <iterator>
#include <algorithm>
#include <functional>
#include <math.h>
#include <random>
#define MIN(a,b) ( (a) <= (b) ? (a) : (b))
#define TEMP 1
#define EPS 0.0001

typedef struct {
    double x,y,z;
} XYZ;


//prototype for function to write Monomer data to file





std::mt19937 generator((unsigned)time(NULL));

template<typename T1>
std::ostream& operator <<( std::ostream& out, const std::vector<T1>& object );

class CSV_out {
    std::string filename;
public:
    CSV_out(std::string filename);
    template <typename T>
    void write_data(T& data, int& sizeN, std::vector<int>& sizeM,std::vector<int>& sizeP, bool O, int samples=100);
};


class Monomer_list {
public:
    double x,y,z;
    int charge, P_number, M_number;
    int ID;
    
    Monomer_list();
    void update(double x1, double y1, double z1, int charge1, int P, int M);
    void updatexyz(double x1, double y1, double z1);
};

//typedefs
typedef boost::multi_array<Monomer_list, 3> MA_datatype;
typedef boost::multi_array<int, 3> LT_datatype;
typedef boost::multi_array_types::index_range range;
typedef std::vector<std::vector<Monomer_list>> VA_datatype;

class MC_mover {
public:
    LT_datatype Clist;
    
    
    VA_datatype Tlist;
    VA_datatype tTlist;
    
    
    std::vector<std::vector<std::vector<int>>> int_vec;
    std::unordered_map<int, std::vector<std::vector<Monomer_list *>>> out_dic;
    std::unordered_map<int, std::vector<std::vector<Monomer_list *>>> tout_dic;
    
    double i_eng;
    double o_eng;
    double t_eng;
    double t_dis;
    int sizeN;
    std::vector<int> sizeP;
    std::vector<int> sizeM;
    int esi=1;
    int change;
    
    MC_mover(std::vector<int> sizeM, std::vector<int> sizeP, int N);
    double hamil(double,int,int);
    double dist(Monomer_list, Monomer_list);
    void update_ints();
    VA_datatype get_Alist();
    LT_datatype get_Clist();
    double tint_eng(VA_datatype);
    double tout_eng(VA_datatype);
    double rint_eng(VA_datatype&,int);
    double rout_eng(int,std::unordered_map<int, std::vector<std::vector<Monomer_list *>>>&);
    void total_dis();
    void mover();
    void create_combi(std::vector<std::vector<int>>& a_vec, int n, int m);
    void create_combi_uti(std::vector<std::vector<int>>& a_vec, std::vector<int>& t_vec, int n,int m, int inter);
    void create_out_dic();
    double test_eng();
    double distp(Monomer_list *, Monomer_list *);
    void create_combi_multi(std::vector<std::vector<std::vector<int>>>& b_vec, std::vector<int> sizeM, int m);
    
};



template <typename T = double>
T get_user_input(std::string a);

void int_grid(LT_datatype& new_LT, int& sizeN);

template <typename T = int>
std::vector<std::vector<int>> ask_charge(std::vector<T>& sizeM,std::vector<T>& sizeP,std::vector<std::vector<T>>& chg_p_monomer);

void initialize_boardh(int& sizeN, std::vector<int>& sizeP, std::vector<int>& sizeM, VA_datatype& M, LT_datatype& L, VA_datatype& M1);

template <typename T = int>
void create_1d_x_yh(T const& sizeN ,std::vector<T> & X,std::vector<T> & Y);

template <typename T>
void get_new_x_yh(T& x, T& y,T& indexs,T& counter, std::vector<T> X, std::vector<T> Y);

void create_combi_multi(std::vector<std::vector<int>>& a_vec, int n, int m, int P);

void create_VA(std::vector<int>& sizeM, std::vector<int>& sizeP,VA_datatype& M, VA_datatype& M1){
    for (int i=0; i!=sizeM.size(); ++i) {
        int k =0;
        while (k<sizeP[i]) {
            k+=1;
            std::vector<Monomer_list> child;
            std::vector<Monomer_list> child0;
            for (int j=0; j!=sizeM[i]; ++j) {
                Monomer_list n;
                Monomer_list n1;
                child.push_back(n);
                child0.push_back(n1);
            }
            M.push_back(child);
            M1.push_back(child0);
        }
    }
}

std::vector<int> input_vec(std::string s){
    std::string line;
    int number;
    std::vector<int> numbers;
    
    std::cout << "Enter "<< s <<" separated by spaces: ";
    std::getline(std::cin, line);
    std::istringstream stream(line);
    while (stream >> number)
        numbers.push_back(number);
    return numbers;
}

int main(int argc, const char * argv[]){
    CSV_out new_file("test");
    std::vector<int> sP = input_vec("polymer");
    std::vector<int> sM = input_vec("monomer");
    int N = get_user_input("Size of box: ");
    int R = get_user_input("MC mover (in multiples of 100): ");
    MC_mover test(sM,sP,N);
//    for(auto& i:test.int_vec){
//        std::cout << "[";
//        for(auto& j:i){
//            std::cout << j << "\n";
//        }
//        std::cout << "]";
//    }
    int samples = get_user_input("How many samples to take: ");
    new_file.write_data(test.Tlist, test.sizeN, test.sizeM, test.sizeP,0,samples +1);
   // new_file.write_data(test.Tlist, test.sizeN, test.sizeM, test.sizeP);
 //   std::cout << test.Tlist[2][2].x << "\n";
    int counter =0;
    for (int i=0; i<R; ++i) {
        if (i%(int((1./samples)*R)) ==0){
            std::cout << "Percent complete: " << (float(i)/float(R))*100 << "%" << "\n";
            if(counter<samples){
                new_file.write_data(test.Tlist, test.sizeN, test.sizeM, test.sizeP,1);
                counter += 1;
            }
        }
        test.mover();
    }
    std::cout << "Percent complete: 100%. Gratz." << "\n";
    std::cout << test.change;
   // std::cout << test.Tlist[2][2].x << "\n";

    return 0;
    
}




template <typename T>
void get_new_x_yh(T& x, T& y,T& indexs,T& counter, std::vector<T> X, std::vector<T> Y){
    if(indexs>=X.size()){
        counter+=1;
        indexs=0;
        x=0;
        y=0;
    } else{
        indexs+=( generator() % (5 - 3) ) + 1;; //change this to initilize random lattice
        if (indexs>=X.size()){
            counter+=1;
            indexs=0;
            x=0;
            y=0;
        } else{
            x=X[indexs];
            y=Y[indexs];
        }
        
    }
    
}

template <typename T>
void create_1d_x_yh(T const& sizeN ,std::vector<T>& X,std::vector<T>& Y){
    int counter_x = 0;
    int counter_y =0;
    int counter_size = 0;
    int i = 0;
    while (i<sizeN*sizeN) {
        X[i]=counter_x;
        counter_x+=1;
        Y[i]=counter_y;
        if(counter_size<sizeN-1){
            counter_size+=1;
        } else{
            counter_size=0;
            counter_y+=1;
            counter_x=0;
        }
        i+=1;
    }
    
}

void initialize_boardh(int& sizeN, std::vector<int>& sizeP, std::vector<int>& sizeM, VA_datatype& M, LT_datatype& L, VA_datatype& M1){
    
    std::vector<std::vector<int>> charge_dist;
    ask_charge(sizeM,sizeP,charge_dist);
   // for(auto& i:charge_dist){
      //  std::cout << i;
   // }
    int it = *max_element(sizeM.begin(), sizeM.end());
    int check_n_m = floor(sizeN/it);
    //int total_P = boost::accumulate(sizeP,0);
    std::vector<int> tX(sizeN*sizeN);
    std::vector<int> tY(sizeN*sizeN);
    
    
    create_1d_x_yh(sizeN, tX, tY);
    
    std::cout << tX << tX.size() <<"\n";
    int counter = 0;
    int x = 0;
    int y = 0;
    int z;
    int indexs = 0;
    for (int i =0; i!=M.size(); ++i) {
    /// this is a temp hack to deal with the initilization problem. Fixe is TODO.
        
        int temp = 0;
        for (int j=0;j!=sizeM.size();++j){
            if (sizeM[j]==M[i].size()){
                temp = j;
            }
            
        }
        std::cout << indexs << "," << temp << "\n";
      //// end of hack.
        z = counter*sizeM[temp];
        std::cout << x << "," <<y << "," <<z <<"\n";
        
        for (int j = 0; j!=M[i].size(); ++j) {
            M[i][j].update(x,y,z,charge_dist[i][j],i+1,j+1);
            M1[i][j].update(x,y,z,charge_dist[i][j],i+1,j+1);
            L[x][y][z]=1.0;
            z+=1;
        }
    
        get_new_x_yh(x,y,indexs,counter,tX,tY);
        if (check_n_m<counter) {
            break;
        }
    }
    
}



void int_grid(LT_datatype& new_LT, int& sizeN){
    for (int i = 0; i< sizeN; ++i) {
        for (int j = 0; j < sizeN; ++j) {
            for (int k = 0; k < sizeN; ++k) {
                new_LT[i][j][k]=0;
            }
        }
    }
}


Monomer_list::Monomer_list(){
    x = 0;
    y = 0;
    z = 0;
    charge = 0;
    P_number = 0;
    M_number = 0;
    ID = 0;
}

void Monomer_list::update(double x1, double y1, double z1, int charge1, int P, int M){
    x = x1;
    y = y1;
    z = z1;
    charge = charge1;
    P_number = P;
    M_number = M;
}

void Monomer_list::updatexyz(double x1, double y1, double z1){
    x = x1;
    y = y1;
    z = z1;
}


template <typename T>
T get_user_input(std::string a){
    T user_input;
    std::cout<< a;
    scanf("%lf", &user_input);
    return user_input;
}






template <typename T>
std::vector<std::vector<int>> ask_charge(std::vector<T>& sizeM,std::vector<T>& sizeP,std::vector<std::vector<T>>& chg_p_monomer){
    
    bool checker = true;

    std::string charge_dist;
    for (int i=0; i!=sizeM.size();++i) {
        checker = true;
        charge_dist.clear();
        while (checker == true) {
            std::cout << "Charge distribution of one polymer which has length " << sizeM[i] << " : ";
            std::cin>>charge_dist;
            if (charge_dist.length() == sizeM[i]*2) {
                checker = false;
            }
        }
        
        std::vector<int> child;
        for(int j = 0; j<charge_dist.length();++j){
            switch (charge_dist[j]) {
                case 'P':
                    child.push_back(charge_dist[j+1] - '0');
                    break;
                case 'N':
                    child.push_back(-(charge_dist[j+1] - '0'));
                    break;
                case 'O':
                    child.push_back(-(charge_dist[j+1] - '0'));
                    break;
            }
        }
        for(int k=0; k!=sizeP[i]; ++k){
            chg_p_monomer.push_back(child);
        }
    }
    
    return chg_p_monomer;
}


CSV_out::CSV_out(std::string filename1){
    filename = filename1;
}

template <typename T>
void CSV_out::write_data(T& coords,int& sizeN, std::vector<int>& sizeM, std::vector<int>& sizeP, bool O, int samples){
    std::fstream file;
    if(O==1){
        file.open(filename, std::ios::out | std::ios::ate | std::ios::app);
    } else{
        file.open(filename, std::ios::out | std::ios::trunc);
    }
    std::ostringstream tstring;
    if(O==0){
        tstring << sizeN << ";";
        
        for(int i=0; i!=sizeP.size(); ++i){
            tstring << sizeP[i];
            if (i!=sizeP.size()-1){
                tstring<< ",";
            }
        }
        tstring << ";";
        for(int i=0; i!=sizeM.size(); ++i){
            tstring << sizeM[i];
            if (i!=sizeM.size()-1){
                tstring<< ",";
            }
        }
        tstring << ";";
        tstring << samples;
        file << tstring.str();
        file << "\n";
    }
    
    for (auto& i : coords) {
        for (auto& j : i) {
            file << j.x <<  "," << j.y << "," << j.z << "," << j.charge;
            file << "\n";
        }
    }
    
    
}

MC_mover::MC_mover(std::vector<int> sM, std::vector<int> sP, int N){
    sizeM=sM;
    sizeP=sP;
    sizeN=N;
    create_VA(sizeM, sizeP, Tlist, tTlist);
    Clist.resize(boost::extents[sizeN][sizeN][sizeN]);
    initialize_boardh(sizeN, sizeP, sizeM, Tlist, Clist, tTlist);

    create_combi_multi(int_vec, sizeM, 2);
    
    create_out_dic();
    i_eng = tint_eng(Tlist);
    o_eng = tout_eng(Tlist);
    t_dis = 0.0;
    t_eng = i_eng + o_eng;
}

double MC_mover::dist(Monomer_list p1, Monomer_list p2){
    double xt = MIN(std::abs(p1.x-p2.x),std::abs(p1.x-p2.x -sizeN));
    double yt = MIN(std::abs(p1.y-p2.y),std::abs(p1.y-p2.y -sizeN));
    double zt = MIN(std::abs(p1.z-p2.z),std::abs(p1.z-p2.z -sizeN));
    double xm = MIN(xt,std::abs(p1.x-p2.x + sizeN));
    double ym = MIN(yt,std::abs(p1.y-p2.y + sizeN));
    double zm = MIN(zt,std::abs(p1.z-p2.z + sizeN));
    double distance = sqrt(pow(xm, 2) + pow(ym, 2) + pow(zm, 2));
    
    return distance;
}

/*
double MC_mover::dist(Monomer_list p1, Monomer_list p2){
    double xm  = std::abs(p1.x-p2.x);
    double xm1 = std::abs(p1.x-p2.x -sizeN);
    double ym  = std::abs(p1.y-p2.y);
    double ym1 = std::abs(p1.y-p2.y -sizeN);
    double zm  = std::abs(p1.z-p2.z);
    double zm1 = std::abs(p1.z-p2.z -sizeN);
    
    double xm11 = std::abs(p1.x-p2.x +sizeN);
    double ym11 = std::abs(p1.y-p2.y +sizeN);
    double zm11 = std::abs(p1.z-p2.z +sizeN);
    
    
    double distance = sqrt(pow(xm, 2) + pow(ym, 2) + pow(zm, 2));
    double distance1 = sqrt(pow(xm1, 2) + pow(ym1, 2) + pow(zm1, 2));
    double distance11 = sqrt(pow(xm11, 2) + pow(ym11, 2) + pow(zm11, 2));
    return MIN(distance,MIN(distance1,distance11));
}
*/
double MC_mover::hamil(double d, int c1, int c2){
    
    double temp=((1./d)*c1*c2 - esi*2./(pow(d,6)) + esi*1./(pow(d,12)) + 0.0*pow(d,4));
    
    return temp;
}

double MC_mover::tint_eng(VA_datatype M){
    double t_e=0.0;
    for (int j = 0; j!=M.size();++j){
        int temp = 0;
        for (int i=0;i!=sizeM.size();++i){
            if (sizeM[i]==M[j].size()){
                temp = i;
            }
            
        }
        
        for (auto& inter:int_vec[temp]){
            t_e += hamil(dist(M[j][inter[0]],M[j][inter[1]]), M[j][inter[0]].charge, M[j][inter[1]].charge);
        }
    }
    return t_e;
}



void MC_mover::create_out_dic(){
    int u_counter=0;
    for (int i = 0; i!=Tlist.size(); ++i) {
        for (int j = 0; j !=Tlist[i].size(); ++j) {
            Tlist[i][j].ID = u_counter;
            tTlist[i][j].ID = u_counter;
            int p =i;
            int counter = 0;
            std::vector<std::vector<Monomer_list *>> pVal;
            std::vector<std::vector<Monomer_list *>> ttpVal;
            while (counter < Tlist.size()-1) {
                int m = 0;
                while (m<Tlist[(p+1)%(Tlist.size())].size()) {
                    Monomer_list *pM1;
                    pM1 = &(Tlist[i][j]);
                    Monomer_list *pM2;
                    pM2= &(Tlist[(p+1)%(Tlist.size())][m]);
                    std::vector<Monomer_list *> tpVal = {pM1,pM2};
                    
                    Monomer_list *tpM1;
                    tpM1 = &(tTlist[i][j]);
                    Monomer_list *tpM2;
                    tpM2 = &(tTlist[(p+1)%(Tlist.size())][m]);
                    std::vector<Monomer_list *> tttpVal = {tpM1,tpM2};
                    
                    pVal.push_back(tpVal);
                    ttpVal.push_back(tttpVal);
                    m+=1;
                }
                p+=1;
                counter += 1;
            }
        //    std::cout << "size of vec: " << pVal.size() << "\n";
            out_dic.emplace(u_counter, pVal);
            tout_dic.emplace(u_counter, ttpVal);
            u_counter+=1;
        }
    }
}


double MC_mover::tout_eng(VA_datatype M){
    double t_e = 0.0;
    for (int i = 0 ; i!=M.size(); ++i) {
        for (int j = 0; j!=M[i].size(); ++j) {
            int p =i;
            // std::vector<std::vector<Monomer_list>> pVal;
            //  std::vector<std::vector<Monomer_list>> ttpVal;
            while (p<M.size()-1) {
                int m = 0;
                while (m<M[i].size()) {
                    // std::cout << "hi" << "\n";
                    //  Monomer_list pM1;
                    //  pM1 = Tlist[i][j][0];
                    //  Monomer_list pM2;
                    //  pM2= Tlist[p+1][m][0];
                    // std::vector<Monomer_list> tpVal = {pM1,pM2};
                    
                    //  Monomer_list tpM1;
                    //  tpM1 = tTlist[i][j][0];
                    //  Monomer_list tpM2;
                    //  tpM2 = tTlist[p+1][m][0];
                    //  std::vector<Monomer_list> tttpVal = {tpM1,tpM2};
                    
                    //  pVal.push_back(tpVal);
                    //  ttpVal.push_back(tttpVal);
                    t_e+=hamil(dist(M[i][j], M[p+1][m]), M[i][j].charge, M[p+1][m].charge);
                    m+=1;
                    
                }
                // std::cout << "," << "\n";
                p+=1;
            }
            // std::cout << ";" << "\n";
            // out_dic.emplace(i*sizeM + j, pVal);
            // tout_dic.emplace(i*sizeM + j, ttpVal);
            //std::cout << "place: " << i*sizeM + j << "\n";
        }
    }
    return t_e;
}





void MC_mover::mover(){
    
    for (int i =0; i!=Tlist.size(); ++i) {
        tTlist = Tlist;
        int tc=( generator() % (Tlist[i].size() - 0) ) + 0;
        int mc=( generator() % (6 - 1 + 1) ) + 1;
        //float r3 = -360 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(360+360)));
        float r = 0.5*static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        //std::cout << r << "\n";
        double& tx= tTlist[i][tc].x;
        double& ty= tTlist[i][tc].y;
        double& tz= tTlist[i][tc].z;
        //std::cout << tx << "," <<ty << "," <<tz << "\n";
        //std::cout << Tlist.size() << "," << i << "\n";
        switch (mc) {
            case 1:
                tx = std::abs(fmod((tx + 1*r + sizeN),sizeN));
                ty = std::abs(fmod((ty + 1*r + sizeN),sizeN));
                tz = std::abs(fmod((tz + 1*r + sizeN),sizeN));
          //      std::cout << tx << "," <<ty << "," <<tz << "\n";
                //std::cout << tTlist[i][tc].x << "," <<tTlist[i][tc].y << "," <<tTlist[i][tc].z << "\n";
                break;
                
            case 2:
                tx = std::abs(fmod((tx + 1*r + sizeN),sizeN));
                ty = std::abs(fmod((ty + 1*r + sizeN),sizeN));
                tz = std::abs(fmod((tz - 1*r + sizeN),sizeN));
           //     std::cout << tx << "," <<ty << "," <<tz << "\n";
               // std::cout << tTlist[i][tc].x << "," <<tTlist[i][tc].y << "," <<tTlist[i][tc].z << "\n";
                break;
            case 3:
                tx = std::abs(fmod((tx + 1*r + sizeN),sizeN));
                ty = std::abs(fmod((ty - 1*r + sizeN),sizeN));
                tz = std::abs(fmod((tz - 1*r + sizeN),sizeN));
           //     std::cout << tx << "," <<ty << "," <<tz << "\n";
                //std::cout << tTlist[i][tc].x << "," <<tTlist[i][tc].y << "," <<tTlist[i][tc].z << "\n";
                break;
            case 4:
                tx = std::abs(fmod((tx - 1*r + sizeN),sizeN));
                ty = std::abs(fmod((ty + 1*r + sizeN),sizeN));
                tz = std::abs(fmod((tz + 1*r + sizeN),sizeN));
             //   std::cout << tx << "," <<ty << "," <<tz << "\n";
                //std::cout << tTlist[i][tc].x << "," <<tTlist[i][tc].y << "," <<tTlist[i][tc].z << "\n";
                break;
            case 5:
                tx = std::abs(fmod((tx - 1*r + sizeN),sizeN));
                ty = std::abs(fmod((ty - 1*r + sizeN),sizeN));
                tz = std::abs(fmod((tz + 1*r + sizeN),sizeN));
              //  std::cout << tx << "," <<ty << "," <<tz << "\n";
               // std::cout << tTlist[i][tc].x << "," <<tTlist[i][tc].y << "," <<tTlist[i][tc].z << "\n";
                break;
            case 6:
                tx = std::abs(fmod((tx - 1*r + sizeN),sizeN));
                ty = std::abs(fmod((ty - 1*r + sizeN),sizeN));
                tz = std::abs(fmod((tz - 1*r + sizeN),sizeN));
              //  std::cout << tx << "," <<ty << "," <<tz << "\n";
                //std::cout << tTlist[i][tc].x << "," <<tTlist[i][tc].y << "," <<tTlist[i][tc].z << "\n";
                break;
        }
        
        bool good_dis = false;
        if (Tlist[i].size() > 1){
            if (tc==0) {
                //  std::cout << tTlist[i][tc].x << ","<< tTlist[i][tc].y << ","<< tTlist[i][tc].z << "," << tTlist[i][tc+1].x<< ","<< tTlist[i][tc+1].y<< ","<< tTlist[i][tc+1].z << "\n";
                // std::cout << dist(tTlist[i][tc],tTlist[i][tc+1]) << "\n";
                if (dist(tTlist[i][tc],tTlist[i][tc+1]) < 1.3) {
                    
                    good_dis = true;
                }
            } else if (tc==Tlist[i].size()-1){
                //    std::cout << tTlist[i][tc].x << ","<< tTlist[i][tc].y << ","<< tTlist[i][tc].z << "," << tTlist[i][tc+1].x<< ","<< tTlist[i][tc+1].y<< ","<< tTlist[i][tc+1].z << "\n";
                
                //  std::cout << dist(tTlist[i][tc],tTlist[i][tc+1]) << "," << dist(tTlist[i][tc],tTlist[i][tc-1]) << "\n";
                if (dist(tTlist[i][tc],tTlist[i][tc-1]) < 1.3) {
                    good_dis = true;
                }
            } else {
                //   std::cout << tTlist[i][tc].x << ","<< tTlist[i][tc].y << ","<< tTlist[i][tc].z << "," << tTlist[i][tc-1].x<< ","<< tTlist[i][tc-1].y<< ","<< tTlist[i][tc-1].z << "\n";
                // std::cout << dist(tTlist[i][tc],tTlist[i][tc-1]) << "\n";
                if (dist(tTlist[i][tc],tTlist[i][tc+1]) < 1.3 && dist(tTlist[i][tc],tTlist[i][tc-1]) < 1.3) {
                    
                    good_dis = true;
                }
            }
        } else {
            good_dis=true;
        }
        
        
        if (good_dis==true) {
            bool lattice_check = false; //((Clist[int(tx)][int(ty)][int(tz)])==1.0);
            
            // std::cout << "Coords: "<< tx <<"," <<ty <<","<<tz<<"\n";
            if(lattice_check==false){
                double riold_eng = rint_eng(Tlist, i);
                double roold_eng = rout_eng(Tlist[i][tc].ID, out_dic);
                double rinew_eng = rint_eng(tTlist, i);
                double ronew_eng = rout_eng(tTlist[i][tc].ID, tout_dic);

                //std::cout << i << "," << tc <<"\n";
                if(rinew_eng == 0 && riold_eng ==0){
                    rinew_eng+=1;
                }
                
                
                
                
                
                // proper MHMC scheme; more compact.
                double min_out_energy = MIN(1.,std::exp(-(ronew_eng-roold_eng)/TEMP));
                double min_in_energy = MIN(1., std::exp(-(rinew_eng-riold_eng)/TEMP));
                
                float r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                
                if (r1<min_in_energy && r1<min_out_energy) {
                    Tlist[i][tc].updatexyz(tx, ty, tz);
                    change+=1;
                    Clist[int(tx)][int(ty)][int(tz)]=1.0;
                }
            
            
                
                
                
                
                
                
                
                
//               // std::cout << riold_eng << "," << roold_eng << "," << rinew_eng << "," << ronew_eng << "\n";
//                if(rinew_eng < riold_eng && ronew_eng < roold_eng){
//                    //std::cout << riold_eng << "," << roold_eng << "," << rinew_eng << "," << ronew_eng << "\n";
//                    Tlist[i][tc].updatexyz(tx, ty, tz);
//                    change+=1;
//                    (Clist[int(tx)][int(ty)][int(tz)])=1.0;
//                } else {
//                    float r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
//                    if (std::exp(-(rinew_eng-riold_eng)/TEMP) > r1 && std::exp(-(ronew_eng-roold_eng)/TEMP) > r1) {
//                        Tlist[i][tc].updatexyz(tx, ty, tz);
//                        change+=1;
//                        Clist[int(tx)][int(ty)][int(tz)]=1.0;
//                    }
                
              //  std::cout <<change << "\n";
            }

        }
        
        
        
        //
    }
    
}


void MC_mover::create_combi_multi(std::vector<std::vector<std::vector<int>>>& b_vec, std::vector<int> sizeM, int m){
    for (auto& i:sizeM) {
        std::vector<std::vector<int>> temp;
        create_combi(temp, i, 2);
        b_vec.push_back(temp);
    }
}


void MC_mover::create_combi(std::vector<std::vector<int>>& a_vec, int n, int m){
    std::vector<int> t_vec;
    create_combi_uti(a_vec, t_vec, n, m, 1);
}


void MC_mover::create_combi_uti(std::vector<std::vector<int>>& a_vec, std::vector<int>& t_vec, int n, int m, int inter){
    
    
    if (m==0) {
        a_vec.push_back(t_vec);
        return;
    }
    
    for (int i = inter; i<=n; ++i) {
        t_vec.push_back(i-1); // changed to allow the combinations to represent indexs starting at 0 and ending at n-1. With n = sizeM.
        create_combi_uti(a_vec, t_vec, n, m-1, i+1);
        t_vec.pop_back();
    }
}




double MC_mover::rint_eng(VA_datatype& M,int i){
    double t_e=0.0;
    int temp = 0;
    for (int j=0;j!=sizeM.size();++j){
        if (sizeM[j]==M[i].size()){
            temp = j;
        }
    
    }
    
    for (auto& inter:int_vec[temp]){
        t_e += hamil(dist(M[i][inter[0]],M[i][inter[1]]), M[i][inter[0]].charge, M[i][inter[1]].charge);
    }
    return t_e;
}

double MC_mover::rout_eng(int i,std::unordered_map<int, std::vector<std::vector<Monomer_list *>>>& dic){
    // int i is the key value of the momomer. Namely, Polymer_number*monomer_count + monomer_number. Where number implies the target particle and count is = sizeM.
    double t_e=0.0;
    std::vector<std::vector<Monomer_list *>> interaction_M = dic[i];
    for (auto & vec : interaction_M) {
     //   std::cout << vec[0]->x << ","<< vec[0]->y << "," <<vec[0]->z << "\n";
        t_e+= hamil(distp(vec[0], vec[1]), vec[0]->charge, (vec[1])->charge);
    }
    return t_e;
}



/*
double MC_mover::distp(Monomer_list * p1, Monomer_list * p2){
    double xm  = std::abs(p1->x-p2->x);
    double xm1 = std::abs(p1->x-p2->x -sizeN);
    double ym  = std::abs(p1->y-p2->y);
    double ym1 = std::abs(p1->y-p2->y -sizeN);
    double zm  = std::abs(p1->z-p2->z);
    double zm1 = std::abs(p1->z-p2->z -sizeN);
    
    double xm11 = std::abs(p1->x-p2->x +sizeN);
    double ym11 = std::abs(p1->y-p2->y +sizeN);
    double zm11 = std::abs(p1->z-p2->z +sizeN);
    
    double distance = sqrt(pow(xm, 2) + pow(ym, 2) + pow(zm, 2));
    double distance1 = sqrt(pow(xm1, 2) + pow(ym1, 2) + pow(zm1, 2));
    double distance11 = sqrt(pow(xm11, 2) + pow(ym11, 2) + pow(zm11, 2));
    return MIN(distance,MIN(distance1,distance11));
}
*/

double MC_mover::distp(Monomer_list * p1, Monomer_list * p2){
    double xt = MIN(std::abs(p1->x-p2->x),std::abs(p1->x-p2->x -sizeN));
    double yt = MIN(std::abs(p1->y-p2->y),std::abs(p1->y-p2->y -sizeN));
    double zt = MIN(std::abs(p1->z-p2->z),std::abs(p1->z-p2->z -sizeN));
    double xm = MIN(xt,std::abs(p1->x-p2->x + sizeN));
    double ym = MIN(yt,std::abs(p1->y-p2->y + sizeN));
    double zm = MIN(zt,std::abs(p1->z-p2->z + sizeN));
    double distance = sqrt(pow(xm, 2) + pow(ym, 2) + pow(zm, 2));
    
    return distance;
}



double MC_mover::test_eng(){
    double t_e;
    t_e = rout_eng(0,out_dic);
    return t_e;
}




// print a vector
template<typename T1>
std::ostream& operator <<( std::ostream& out, const std::vector<T1>& object )
{
    out << "[";
    if ( !object.empty() )
    {
        for(typename std::vector<T1>::const_iterator
            iter = object.begin();
            iter != --object.end();
            ++iter) {
            out << *iter << ", ";
        }
        out << *--object.end();
    }
    out << "]";
    return out;
}


int LineLineIntersect(
                      XYZ p1,XYZ p2,XYZ p3,XYZ p4,XYZ *pa,XYZ *pb,
                      double *mua, double *mub)
{
    
    XYZ p43;
    p43.x = p4.x - p3.x;
    p43.y = p4.y - p3.y;
    p43.z = p4.z - p3.z;
    
    if (std::abs(p43.x) < EPS && std::abs(p43.y) < EPS && std::abs(p43.z) < EPS){
        return(false);
    }
    XYZ p13;
    p13.x = p1.x - p3.x;
    p13.y = p1.y - p3.y;
    p13.z = p1.z - p3.z;
    
    XYZ p21;
    p21.x = p2.x - p1.x;
    p21.y = p2.y - p1.y;
    p21.z = p2.z - p1.z;
    
    if (std::abs(p21.x) < EPS && std::abs(p21.y) < EPS && std::abs(p21.z) < EPS){
        return(false);
    }
    double d1343,d4321,d1321,d4343,d2121;
    d1343 = p13.x * p43.x + p13.y * p43.y + p13.z * p43.z;
    d4321 = p43.x * p21.x + p43.y * p21.y + p43.z * p21.z;
    d1321 = p13.x * p21.x + p13.y * p21.y + p13.z * p21.z;
    d4343 = p43.x * p43.x + p43.y * p43.y + p43.z * p43.z;
    d2121 = p21.x * p21.x + p21.y * p21.y + p21.z * p21.z;
    
    double denom;
    denom = d2121 * d4343 - d4321 * d4321;
    if (std::abs(denom) < EPS)
    {
        return(false);
    }
    
    double numer;
    numer = d1343 * d4321 - d1321 * d4343;
    
    *mua = numer / denom;
    *mub = (d1343 + d4321 * (*mua)) / d4343;
    
    pa->x = p1.x + *mua * p21.x;
    pa->y = p1.y + *mua * p21.y;
    pa->z = p1.z + *mua * p21.z;
    pb->x = p3.x + *mub * p43.x;
    pb->y = p3.y + *mub * p43.y;
    pb->z = p3.z + *mub * p43.z;
    
    return(true);
}


