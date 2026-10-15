#include<iostream>
#include<string.h>
#include<fstream>
#include <unistd.h> 

//int main()
int main(int argc, char *argv[])		// fix by BM
{
  char file[5], string[500];
  std::ifstream in_file;
  std::ofstream out_file;
  int i, xmin, xmax, ymin, ymax;

  int c;
  char *filename = 0;
  /* parse command line arguments using getopt - Added by BM */
  // Syntatx is bhead -i inputfile
  while ( (c = getopt(argc, argv, "i:")) != -1) {
        switch (c) {
        case 'i':
            //printf ("input file '%s'\n", optarg);
            filename = optarg;
            break;
        default:
            printf ("?? getopt returned character code 0%o ??\n", c);
        }
  }
  if (optind < argc) {
        printf ("Invalid argument: ");
        while (optind < argc)
            printf ("%s ", argv[optind++]);
        printf ("\n");
  }


  //in_file.open("b.ps");
  in_file.open(filename);
  in_file >> string >> xmin >> string;
  in_file >> string >> xmax >> string;
  in_file >> string >> ymin >> string;
  in_file >> string >> ymax >> string;
  in_file.close();
  out_file.open("bbox.ps");
  out_file << "%!PS-Adobe-2.0 EPSF-2.0" << '\n';
  out_file << "%%BoundingBox: " << xmin << " ";
  out_file << ymin << " " ;
  out_file << xmax << " " ;
  out_file << ymax << '\n';
  out_file.close();
  return(0);
}
