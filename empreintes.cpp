

#include <fstream>
#include <iostream>
#include <string>
#include <GraphMol/GraphMol.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/stream.hpp>

#include <ctime>
using namespace std;
#define NB_MOLECULES 90130
#define ECHANTILLON 10000

int main( int argc , char **argv ) {


  RDKit::ROMOL_SPTR mol;
  std::string sdf_file = "rdkit_chebi.sdf";
  bool takeOwnership = true;
  RDKit::SDMolSupplier mol_supplier( sdf_file , takeOwnership );
  RDKit::ROMOL_SPTR tabChebi[int(mol_supplier.length())];
  RDKIT_FINGERPRINTS_EXPORT ExplicitBitVect *tabChebiFingerprint[NB_MOLECULES];
  int tableau_molecule[ECHANTILLON];

  float temps;
  clock_t t1, t2;

  string const fichier1 { "time.data" };
  ofstream monFlux { fichier1 };
  if(monFlux)
  {
    for( int i = 0;i < int(mol_supplier.length());i++)
    {
      cout << i << endl;
      RDKit::ROMOL_SPTR mol( mol_supplier[i] );
      if(mol)
      {
        //fingerprint + time
          tabChebi[i] = mol;
          t1 = clock();
          tabChebiFingerprint[i] = RDKFingerprintMol(*tabChebi[i]);//daylight fingerprint
          t2 = clock();
          temps = (float)(t2-t1)/CLOCKS_PER_SEC;
          monFlux << mol->getProp<std::string>( "_Name" ) << " "<< temps << endl;
      }
    }
  }
  else
  {
    cout << "ERREUR: Impossible d'ouvrir le fichier." << endl;
  }
  monFlux.close();

  string const fichier2 { "indice_molecules.data" };
  ifstream monFlux2{ fichier2};
  int molecule,i,j;
  if(monFlux2)
  {


    for( i = 0; i < ECHANTILLON;i++)
    {
      monFlux2 >> molecule;
      tableau_molecule[i] = molecule;
    }
    monFlux2.close();
  }
  else
  {
      cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
  }



  //similarity Tanimoto
  double sim;
  string const fichier3 { "similarity_tanimoto.data" };
  ofstream monFlux3 { fichier3 };
  string const fichier4 { "timec_tanimoto.data" };
  ofstream monFlux4 { fichier4 };
  if(monFlux3 && monFlux4)
  {
    for( i = 1; i < ECHANTILLON ;i++)
    {
      cout <<i<< endl;
      if( tabChebi[tableau_molecule[i]] )
      {
        for( j = 0; j < i ;j++)
        {
          if(tabChebi[tableau_molecule[j]] )
          {
            //std::cout <<  tabChebi[tableau_molecule[i]]->getProp<std::string>( "_Name" ) << "sim "<<  tabChebi[tableau_molecule[j]]->getProp<std::string>( "_Name" )<< std::endl;
            t1 = clock();
            sim = CalcBitmapTanimoto((const unsigned char*)tabChebiFingerprint[tableau_molecule[i]], (const unsigned char*)tabChebiFingerprint[tableau_molecule[j]],tabChebiFingerprint[tableau_molecule[j]]->getNumBits() );
            t2 = clock();
            temps = (float)(t2-t1)/CLOCKS_PER_SEC;
            monFlux3 << sim << " " ;
            monFlux4 << temps << " " ;
            //std::cout <<  " sim =  " << sim  << std::endl;
          }
        }
      }
      monFlux3 << endl;
      monFlux4 << endl;
    }
    monFlux3.close();
    monFlux4.close();
  }
  else
  {
    cout << "ERREUR: Impossible d'ouvrir le fichier en ecriture." << endl;
  }


  return 0;
}
