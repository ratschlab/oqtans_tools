// Authors: Bettina Hepp, Uta Schulze, Cheng Soon Ong, Fabio De Bona, Gunnar Raetsch, Geraldine Jean, Soeren Sonnenburg
// Copyright (C) 2005-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

#include "fill_matrix.h"
using namespace std;

void print_align(Pre_score* matrix, int length_est,  int length_dna,Align_pair* vektor, int result_length, int print_matrix)
{
  //Druck der matrix
  switch(print_matrix)
    {
    case 1:
      {
	cout <<"Matrix: \n";
	
	for (int k = 0; k<length_est; k++)//k Zeile, l Spalte
	  {
	    cout << "Zeile " << k << ":\t"; 
	    for (int l= 0; l<length_dna; l++)
	      {
		cout << ((Pre_score*)matrix + (k*length_dna) + l)->value <<"\t";  
	      } 
	    cout <<"\n";
	  }
	break;	
      }
      
    case 2:
      {
	cout <<"Matrix: \n";
	
	for (int k = 0; k<length_est; k++)//k Zeile, l Spalte
	  {
	    cout << "Zeile " << k << ":\t"; 
	    for (int l= 0; l<length_dna; l++)
	      {
		if (((((Pre_score*)matrix + (k*length_dna) + l)->prev_i) == k-1) && ((((Pre_score*)matrix + (k*length_dna) + l)->prev_j) == l-1))
		  cout << "\\" <<((Pre_score*)matrix + (k*length_dna) + l)->value <<"\t";  
		else if (((((Pre_score*)matrix + (k*length_dna) + l)->prev_i) == k) && ((((Pre_score*)matrix + (k*length_dna) + l)->prev_j) == l-1))
		  cout << "<"<< ((Pre_score*)matrix + (k*length_dna) + l)->value <<"\t";  
		else if (((((Pre_score*)matrix + (k*length_dna) + l)->prev_i) == k -1) && ((((Pre_score*)matrix + (k*length_dna) + l)->prev_j) == l))
		  cout << "|"<< ((Pre_score*)matrix + (k*length_dna) + l)->value <<"\t";  
		else 	   	 { cout << ((Pre_score*)matrix + (k*length_dna) + l)->value << "("; 
		cout << ((Pre_score*)matrix + (k*length_dna) + l)->prev_i << "," << ((Pre_score*)matrix + (k*length_dna) + l)->prev_j << ")\t";  
		}
		
	      } 
	    cout <<"\n";
	  }
        
	cout << "length_est" << length_est << "\n";
        cout << "length_dna" << length_dna << "\n";
        cout << "result_length: " << result_length << "\n";
	
	
	break;	
      } 
      
    case 3:
      {
	//Druck des Alignment
	
	// cout << "Alignment\n" ;
	int j = 0;
	// cout << "resultlength: " << result_length <<"\n";
	while (j<result_length) {
	  int pacman = j; 	
	  cout << "D:(" <<j<<")\t";
	  for (int i=0; ((j<result_length) && (i<50)); i++)
	    {
	      cout << vektor[j].dna_char << " ";	
	      j++;	
	    }
	  
	  cout << "\n";	
	  
	  cout << "E:(" <<pacman<<")\t";
	  for (int i=0;((pacman<result_length) && (i<50)); i++)
	    {
	      cout << vektor[pacman].est_char << " "; 
	      pacman++;
	    }
	  cout << "\n\n";	
	}		
	
	//Ende Druck Alignment
	break ;
      }
      
    default:
      {
	cout <<"";
      }
    }
  
  
  
}

