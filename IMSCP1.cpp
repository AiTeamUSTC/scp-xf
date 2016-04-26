#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<limits.h>
#include<float.h>
#include<memory.h>
#include<string.h>
//#include<sys/times.h>
#include<time.h>
#include<unistd.h>
#include<memory.h>
#include<set>
#define MAXC 50000
#define MAXR 500000
#include<fstream> //利用c++文件读写操作 
#include<iostream>
#include<algorithm>
using namespace std;
typedef struct column{
	int config;
	int score;   
	int time_stamp;
	int is_in_c;
	int cost;
	vector<int> vec_rs;// rows it can cover
}column;

typedef struct row{
	int weight;
	int index;
	vector<int> vec_cs; // cols have can cover this row
}row;
/*声明全局变量*/
int c_num;
int r_num;
column cs[MAXC];
row rs[MAXR];
double alpha[MAXC];

vector<int> uncov_r;
vector<int> best_sol;
int best_solution[200];
set<int> tabu_list;
int tabu[20];
int c_size;
int seed;
int BEST;
int step;
int total_step;
int best_array[MAXC];
int Msize=0;
int migrate_size;
char outname[30];
char recordname[30];
char solname[30];
int island_num;
int BSK=1000;
int BSK_size;
 /*声明函数*/
void init();
void init_best(char* filename);
void localsearch(char* filname);
void record(char* filename);
void add(int);
void remove(int);
int find_best_in_c(int);
void uncov_r_weight_inc();
void update_best_sol();
void migrate_best_sol();
int fitness();
int best_value;
void build_instance2(char *file);
int check();
void print_sol(char* filename);
/*主函数*/
int main(int argc, char *argv[])
{

  
  
  
  for(int i=1;i<=10;i++)
  {
 
    char filename[]="scp41.txt"; 
    best_sol.clear();
    uncov_r.clear();
    cout<<"Please choose island: "<<endl;
    cin>>island_num;
    build_instance2(filename);
    BEST=1;
    switch(island_num)
    {
    	 case 1:
		 {
		    seed=10; strcpy(outname,"island1.txt");  strcpy(recordname,"record1.txt");  strcpy(solname,"soultion1.txt");break;
		 }
	     case 2:
		 {
		 	seed=11; strcpy(outname,"island2.txt");  strcpy(recordname,"record2.txt");  strcpy(solname,"soultion2.txt");break;
		 }
	     case 3:
		 {
		 	seed=12; strcpy(outname,"island3.txt");  strcpy(recordname,"record3.txt");  strcpy(solname,"soultion3.txt");break;
		 }
	     case 4:
		 {
		 	seed=13; strcpy(outname,"island4.txt");  strcpy(recordname,"record4.txt");  strcpy(solname,"soultion4.txt");break;
		 }
         case 5:
         {
         	seed=14; strcpy(outname,"island5.txt");  strcpy(recordname,"record5.txt");  strcpy(solname,"soultion5.txt");break;
		 }
		 case 6:
		 {
		    seed=15; strcpy(outname,"island6.txt");  strcpy(recordname,"record6.txt");  strcpy(solname,"soultion6.txt");break;
		 }
	     case 7:
		 {
		 	seed=16; strcpy(outname,"island7.txt");  strcpy(recordname,"record7.txt");  strcpy(solname,"soultion7.txt");break;
		 }
	     case 8:
		 {
		 	seed=17; strcpy(outname,"island8.txt");  strcpy(recordname,"record8.txt");  strcpy(solname,"soultion8.txt");break;
		 }
	     case 9:
		 {
		 	seed=18; strcpy(outname,"island9.txt");  strcpy(recordname,"record9.txt");  strcpy(solname,"soultion9.txt");break;
		 }
         case 10:
         {
         	seed=19; strcpy(outname,"island10.txt");  strcpy(recordname,"record10.txt");  strcpy(solname,"soultion10.txt");break;
		 }
        
		
         
    }
    total_step=100000;
    srand(seed);
    ofstream out;
	init(); 
//	times(&start);
    clock_t start,end;
	start=clock(); 
	init_best(outname);
	int mode; 
	cout<<"Please choose MODE 1:Record  MODE 2:Search: "<<endl;
	cin>>mode;
	switch(mode)
	{  
	   case 1:
	   {
	   	 record(outname);break;
	   }
	
       case 2:
       {
		localsearch(outname);break;
	   }
    }
//	times(&end);
    end=clock();
    out.open(outname,ios::app);
	printf("final best:%d\n",best_value);
	out<<"final best:"<<best_value<<endl;
	printf("search steps is %d, seed is %d\n",step,seed);
	out<<"search steps is "<<step<<","<<"seed is "<<seed<<endl;
	if(!check()) 
	{
		printf("wrong answer \n");
		out<<"wrong anwser"<<endl;
	}
	printf("time cost %f s\n", (end-start)*1.0/CLOCKS_PER_SEC);
	out<<"time cost "<<(double)((end-start)*1.0/CLOCKS_PER_SEC)<<endl;
	out.close();
	/*switch(mode)
	{
		case 1:
		{
			print_sol(recordname);
			break;
		} 
		case 2: 
		{
			print_sol(solname);
			break;
		}
	}*/
	if (best_value<BSK)
    {
    	BSK=best_value;
    	BSK_size=best_sol.size();
    	memset(best_solution,0,sizeof(best_solution));
    	for(int i=0;i<best_sol.size();i++)
    	{
    		best_solution[i]=best_sol[i];                              
		}
	}   
 }
	/*int B[1000];
	for (int i=1;i<=r_num;i++)
	{
		B[i]=rs[i].weight;
	}
	ofstream sol;
	sol.open("RW_seed10.txt");
	for(int i=1;i<=r_num;i++)
	{
	  	sol<<"R"<<i<<": "<<B[i]<<endl;                                                                 
	} 
	sol.close();*/
    ofstream bskfile;
    bskfile.open("scp41.bsk.txt");
    bskfile<<"best solution: "<<BSK<<endl;
    for(int i=0;i<BSK_size;i++)
    {
    	bskfile<<best_solution[i]<<endl;
	}
	bskfile.close();
	return 0;
}

void build_instance2(char *file){
	int i,j,k,t;
		
	ifstream f1;
	f1.open(file);
	f1>>r_num>>c_num;
	for(i=1;i<=c_num;i++){
		f1>>k;
		cs[i].cost=k;
		cs[i].vec_rs.clear();
	}
	for(i=1;i<=r_num;i++)
	{
		rs[i].vec_cs.clear();
	}
	for(i=1;i<=r_num;i++)
	{
		f1>>k;
		rs[i].weight=1;
		for(j=0;j<k;j++)
		{
			f1>>t;
			rs[i].vec_cs.push_back(t);
			cs[t].vec_rs.push_back(i);
		}
		uncov_r.push_back(i);
	}
	for(i=1;i<=c_num;i++){
		cs[i].score=0;
		cs[i].config=1;
		cs[i].time_stamp=0;
		cs[i].is_in_c=0;
		for(j=0;j<cs[i].vec_rs.size();j++){
			k=cs[i].vec_rs[j];
			cs[i].score+=rs[k].weight;
		}
		alpha[i]=cs[i].score*1.0/cs[i].cost;
	}
	f1.close();
}
void init(){
	best_sol.clear();
	c_size=0;
	tabu[0]=tabu[1]=tabu[2]=0;
	best_value=INT_MAX;
}
void init_best(char* filename){
	int cnt;
	int i,j,k;
	double max;
	/*构造一个初始解*/
	while(!uncov_r.empty()){
		max=-DBL_MAX;
		cnt=0;
		memset(best_array,0,sizeof(best_array));
		/*找出score最大的列*/ 
		for(i=1;i<=c_num;i++){
			if(cs[i].is_in_c) continue;
			if(max<alpha[i]){
				max=alpha[i];
				best_array[0]=i;
				cnt=1;
			} 
			else if(max==alpha[i]){
				best_array[cnt++]=i;
			}
		}
	//	printf("cnt==%d, max==%f, ",cnt, max);
	/*随机抽取分数最大列，略繁琐*/
		if(cnt>0){
			vector<int> bestA;                           
			int minCost=INT_MAX;
			k=0;
			for(j=0;j<cnt;j++)
			{
				if(best_array[j]< minCost)
				{
					minCost=best_array[j];
					bestA.clear();                                                                                       
					bestA.push_back(j);  
				} 
				else if( minCost==best_array[j]) 
				{
					bestA.push_back(j);
				}
			}
			j=rand()%bestA.size();
			add(best_array[bestA[j]]);
		}
	}
	update_best_sol();
	ofstream out;
	out.open(filename);
	printf("initial solution:%d\n",best_value);
	out<<"initial solution:"<<best_value<<endl;
	out.close();
	
}
void add(int c){
	cs[c].is_in_c=1;
	cs[c].score=-cs[c].score;
	alpha[c]=-alpha[c];
	int i,j,k,t,cnt,s;
	for(i=0;i<cs[c].vec_rs.size();i++){
		j=cs[c].vec_rs[i];
		cnt=0;
		for(k=0;k<rs[j].vec_cs.size();k++){
			t=rs[j].vec_cs[k];
			if(t==c) continue;
			if(cs[t].is_in_c){
				s=t;
				cnt++;
			}
			cs[t].config=1;
		}
		if(cnt==0){ // c is the first one can cover this row 
			for(k=0;k<rs[j].vec_cs.size();k++){
				t=rs[j].vec_cs[k];
				if(t==c) continue;
				cs[t].score-=rs[j].weight;
				alpha[t]=cs[t].score*1.0/cs[t].cost;
			}
		} else if(cnt==1){// c is second one can cover this row
			cs[s].score+=rs[j].weight;
			alpha[s]=cs[s].score*1.0/cs[s].cost;
		}
	}
	
	vector<int>::iterator it=uncov_r.begin();
	while(it!=uncov_r.end()){
		for(i=0;i<rs[*it].vec_cs.size();i++){
			if(rs[*it].vec_cs[i]==c){
				break;
			}
		}
		if(i<rs[*it].vec_cs.size()){
			it=uncov_r.erase(it);
		} else {
			it++;
		}
	}

	c_size++;
}
void remove(int c){
	cs[c].is_in_c=0;
	cs[c].score=-cs[c].score;
	alpha[c]=-alpha[c];
	cs[c].config=0;
	int i,j,k,t,cnt,s;
	for(i=0;i<cs[c].vec_rs.size();i++){
		j=cs[c].vec_rs[i];
		cnt=0;
		for(k=0;k<rs[j].vec_cs.size();k++){
			t=rs[j].vec_cs[k];
			if(t==c) continue;
			if(cs[t].is_in_c){
				cnt++;
				s=t;
			}
			cs[t].config=1;
		}
		if(cnt==0){
			uncov_r.push_back(j);
			for(k=0;k<rs[j].vec_cs.size();k++){
				t=rs[j].vec_cs[k];
				if(t==c) continue;
				cs[t].score+=rs[j].weight;
				alpha[t]=cs[t].score*1.0/cs[t].cost;
			}
		} else if(cnt==1){
			cs[s].score-=rs[j].weight;
			alpha[s]=cs[s].score*1.0/cs[s].cost;
		}
	}
	c_size--;
}
/*检测是否在禁忌表中*/ 
int in_tabu(int i){
	int j;
	int tabu_len=1;
	for(j=0;j<tabu_len;j++){
		if(i==tabu[j]) return 1;
	}
	return 0;
}
/*选择C中要remove的一列*/ 
int find_best_in_c(int allowTabu)
{
	int i, maxc, cnt=0;
	double max=-DBL_MAX;
	for(i=1;i<=c_num;i++)
	{
		if(!cs[i].is_in_c) continue;
		if(allowTabu&&in_tabu(i)) continue;
		if(max<alpha[i])
		{
			max=alpha[i];
			maxc=i;
			best_array[0]=i;
			cnt=1;
		} 
		else if(max==alpha[i])
		{
			best_array[cnt++]=i;
			/*分数相同的情况下比较时间戳，先进先出*/ 
			if(cs[maxc].time_stamp>cs[i].time_stamp)
			{
				maxc=i;
			}
		}
	}
	return maxc;
}
/*对于未被覆盖的行的权值进行调整*/ 
void uncov_r_weight_inc(){
	int i,j,k,total=0,cnt,s;
	for(i=0;i<uncov_r.size();i++){
		rs[uncov_r[i]].weight+=1;
		
		/*相应的可以覆盖该行的列的分数增加*/ 
		for(j=0;j<rs[uncov_r[i]].vec_cs.size();j++){
			k=rs[uncov_r[i]].vec_cs[j];
			cs[k].score+=1;
			
			alpha[k]=cs[k].score*1.0/cs[k].cost;
		}
	}
}
void localsearch(char*filename){
	step=1;
	int i,j,k;
	int best_in_c;
	int maxc;
	double max=0;
	int fit;
    migrate_size=2; 
	ofstream out;
	out.open(filename,ios::app);
	Msize=0;
	//ofstream f2;
	//f2.open("fitness1.txt",ios::app);
	while(step<total_step)
	{
	  
	    //fit=fitness();
		//f2<<fit<<endl;  	    
		if(uncov_r.empty())
		{	
		   /*定义迁移标准*/		
		   Msize++; 
		   if(Msize<migrate_size)
		     update_best_sol();
		    else 
		    {
		  	    migrate_best_sol();
		    }
		    //update_best_sol();
			//cout<<fit<<endl; 
			if(best_value==BEST){
				printf("optimum found\n");
                out<<endl<<"optimum found"<<endl;
                out.close();
				break;
			}
			maxc=find_best_in_c(0);
		//	cout<<"12345566: "<<maxc<<endl; 
			remove(maxc);
			continue;
		}
		/*allow Tabu，使用禁忌表*/ 
		best_in_c=find_best_in_c(1);
	    
		remove(best_in_c);
		cs[best_in_c].time_stamp=step;
		
		while(!uncov_r.empty()){
			i=rand()%uncov_r.size();
			j=uncov_r[i];
		//	printf("uncov size:%d j==%d, rs[j]:%d.",uncov_r.size(),j,rs[j].vec_cs.size());
		    max=-DBL_MAX; 
			for(i=0; i<rs[j].vec_cs.size();i++){
				k=rs[j].vec_cs[i];
				if(!cs[k].config) continue;
				if(max<alpha[k]){
					max=alpha[k];
					maxc=k;
				} else if(max==alpha[k]&&cs[maxc].time_stamp>cs[k].time_stamp){
					maxc=k;
				}
			}
			
			if(cs[maxc].cost+fitness()>=best_value) break;
			add(maxc);
			tabu[0]=maxc;
			cs[maxc].time_stamp=step;
			uncov_r_weight_inc();
		}
        
		step++;
	}
	
}
void record(char* filename)
{
	step=1;
	int i,j,k;
	int best_in_c;
	int maxc;
	double max=0;
	int fit;
    migrate_size=2; 
	ofstream out;
	out.open(filename,ios::app);
	Msize=0;
	while(step<total_step)
	{
	    if(Msize==migrate_size)
	       break;
		if(uncov_r.empty())
		{	 
		     Msize++;
		     update_best_sol();
		     
			if(best_value==BEST)
			{
				printf("optimum found\n");
                out<<endl<<"optimum found"<<endl;
                out.close();
				break;
			}
			maxc=find_best_in_c(0);
			remove(maxc);
			continue;
		}
		/*allow Tabu，使用禁忌表*/ 
		best_in_c=find_best_in_c(1);
	    
		remove(best_in_c);
		cs[best_in_c].time_stamp=step;
		
		while(!uncov_r.empty())
		{
			i=rand()%uncov_r.size();
			j=uncov_r[i];
		//	printf("uncov size:%d j==%d, rs[j]:%d.",uncov_r.size(),j,rs[j].vec_cs.size());
		    max=-DBL_MAX; 
			for(i=0; i<rs[j].vec_cs.size();i++)
			{
				k=rs[j].vec_cs[i];
				if(!cs[k].config) continue;
				if(max<alpha[k])
				{
					max=alpha[k];
					maxc=k;
				} else if(max==alpha[k]&&cs[maxc].time_stamp>cs[k].time_stamp)
				{
					maxc=k;
				}
			}
			
			if(cs[maxc].cost+fitness()>=best_value) break;
			add(maxc);
			tabu[0]=maxc;
			cs[maxc].time_stamp=step;
			uncov_r_weight_inc();
		}
        
		step++;
	}
	
}
int fitness(){
	int s=0;
	int i;
	for(i=1;i<=c_num;i++){
		if(cs[i].is_in_c){
			s+=cs[i].cost;
		}
	}
	return s;
}
void update_best_sol(){
	int i;
	int k=fitness();
	//printf("update: %d ",k);
	if(k<best_value){
		best_sol.clear();
		best_value=k;
		for(i=1;i<=c_num;i++){
			if(cs[i].is_in_c){
				best_sol.push_back(i);
			}
		}
	}
}
void migrate_best_sol(){
	 ifstream Msol;
	 int k;
	 int Sumcost;
	 int size;
	 best_sol.clear();
	 switch(island_num)
	 {
	 	case 1:
	 		{
	 		  Msol.open("record10.txt");
			   break;	
			}
		case 2:
	 		{
	 		  Msol.open("record1.txt");
			   break;	
			}
		case 3:
	 		{
	 		  Msol.open("record2.txt");
			   break;	
			}
		case 4:
	 		{
	 		  Msol.open("record3.txt");
			   break;	
			}
		case 5:
	 		{
	 		  Msol.open("record4.txt");
			   break;	
			}
		case 6:
	 		{
	 		  Msol.open("record5.txt");
			   break;	
			}
		case 7:
	 		{
	 		  Msol.open("record6.txt");
			   break;	
			}
		case 8:
	 		{
	 		  Msol.open("record7.txt");
			   break;	
			}
		case 9:
	 		{
	 		  Msol.open("record8.txt");
			   break;	
			}
		case 10:
	 		{
	 		  Msol.open("record9.txt");
			   break;	
			}
			
	 }
	 
	 Msol>>Sumcost;
	 Msol>>size;
	 if(Msize==migrate_size)
	 {
	    cout<<"migrate"<<endl;
	    if(best_value<Sumcost)
	     update_best_sol();
	    else if(best_value>Sumcost)
	    {   cout<<"migrate_1"<<endl;
	        best_sol.clear();
		    best_value=Sumcost;
	    	for(int i=1;i<=c_num;i++)
	        {
	 	      cs[i].is_in_c=0;
	        }
	        for(int i=1;i<=size;i++)
	        {
	            Msol>>k;
	   	        best_sol.push_back(k);
	 	        cs[k].is_in_c=1;
	        }
	        Msol.close();
		}
		else if(best_value==Sumcost)
		{
			int k=rand()%2;
			if(k==1)
			    update_best_sol();
			else if(k==0)
			{   cout<<"migrate_2"<<endl;
			    best_sol.clear();
		        best_value=Sumcost;
				for(int i=1;i<=c_num;i++)
	            {
	 	          cs[i].is_in_c=0;
	            }
	            for(int i=1;i<=size;i++)
	            {
	              Msol>>k;
	   	          best_sol.push_back(k);
	 	          cs[k].is_in_c=1;
	            }
	            Msol.close();
			}
			
		}
	    
	 }
	 else
	    update_best_sol();
	 
}
int check(){ // check if the solution is a correct cover
	int i,j,k;
	int t[MAXR];
	memset(t,0,MAXR*sizeof(int));
	for(i=0;i<best_sol.size();i++){
		j=best_sol[i];
		for(k=0;k<cs[j].vec_rs.size();k++){
			t[cs[j].vec_rs[k]]=1;
		}
	}
	for(i=1;i<=r_num;i++){
		if(!t[i]) return 0;
	}
	return 1;
}
void print_sol(char*filename){
	ofstream out;
	out.open(filename);   
	out<<fitness()<<endl;   
	out<<best_sol.size()<<endl;                                                                                                                                                                                                                          
	for(int i=0;i<best_sol.size();i++){
		printf("%d,",best_sol[i]);
		out<<best_sol[i]<<endl;
		
	}
	printf("\n");
	out<<endl;

}
