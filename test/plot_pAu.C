void plot_pAu()
{

  ifstream fin("full_data.txt");

  string system;
  string junk;
  double v22raw[30];
  double ev22raw[30];
  double v22sub[30];
  double ev22sub[30];
  for ( int i = 0; i < 30; ++i )
    {
      fin>>system>>junk>>junk>>v22raw[i]>>junk>>ev22raw[i]>>junk>>junk>>junk>>v22sub[i]>>junk>>ev22sub[i];
      cout << i << " " << system << " " << v22raw[i] << " " << ev22raw[i] << " " << v22sub[i] << " " << ev22sub[i] << endl;
    }

}
