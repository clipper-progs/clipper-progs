/*! \file buccaneer-tidy.cpp buccaneer library */
/* (C) 2010 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-tidy.h"
#include "buccaneer-prot.h"

#include <algorithm>

#ifdef scr1 // defined on Windows
# undef scr1
# undef scr2
#endif

bool ModelTidy::tidy( clipper::MiniMol& mol, const clipper::MiniMol& mol_mr, const clipper::MMoleculeSequence& seq ) const
{
  std::vector<int> seqnums = chain_renumber( mol, seq );
  std::vector<int> chnnums = chain_assign( mol, mol_mr, seqnums, rmsd_, nmin_ );
  chain_move( mol, mol_mr, chnnums );
  return update_model( mol, mol_mr, chnnums );
}


std::vector<int> ModelTidy::chain_assign( const clipper::MiniMol& mol, const clipper::MiniMol& mol_mr, const std::vector<int> seqnums, const double rmsd, const int nmin )
{
  if ( mol_mr.size() > 0 )
    return assign_chains( mol, mol_mr );
  else
    return assign_chains( mol, seqnums, rmsd, nmin );
}


bool ModelTidy::chain_move( clipper::MiniMol& mol, const clipper::MiniMol& mol_mr, std::vector<int> chnnums )
{
  const clipper::Spacegroup& spgr = mol.spacegroup();
  const clipper::Cell&       cell = mol.cell();

  if ( mol_mr.size() > 0 ) {

    // We have an MR model, shift all the chains onto it
    ProteinTools::symm_match( mol, mol_mr );

  } else {

    // No MR model, group chains around largest fragement for each source
    // current number of chains
    const int nchn = mol.size();
    // ideal number of chains (chnnums)
    int nsrc = 0;
    for ( int c = 0; c < chnnums.size(); c++ )
      nsrc = std::max( nsrc, chnnums[c]+1 );
    // for each source in turn, group around centre of mass
    for ( int s = 0; s < nsrc; s++ ) {
      clipper::MPolymer mptmp;
      for ( int c = 0; c < nchn; c++ )
	if ( chnnums[c] == s ) {
	  if ( mptmp.size() > 0 )
	    move_chain( mol[c], mptmp, spgr, cell );
	  for ( int r = 0; r < mol[c].size(); r++ )
	    mptmp.insert( mol[c][r] );
	}
    }

  }
  return true;
}


std::vector<int> ModelTidy::assign_chains( const clipper::MiniMol& mol, const clipper::MiniMol& mol_mr )
{
  const int nchn = mol.size();
  const int nsrc = mol_mr.size();

  if ( nchn == 0 ) return std::vector<int>();

  // set up chain numbers
  const double rad = 5.0;
  clipper::MAtomNonBond nb( mol_mr, rad );
  std::vector<int> srcnums( nsrc );
  for ( int c2 = 0; c2 < nsrc; c2++ ) srcnums[c2] = c2;

  // count references
  std::vector<int> chnnums( nchn, -1 );
  for ( int c1 = 0; c1 < nchn; c1++ ) {
    std::vector<int> refcount =
      count_contacts( nb, mol_mr, srcnums, mol[c1], rad );
    int c = 0;
    for ( int c2 = 0; c2 < nsrc; c2++ )
      if ( refcount[c2] > refcount[c] ) c = c2;
    if ( refcount[c] > 5 ) chnnums[c1] = c;
  }

  return chnnums;
}


std::vector<int> ModelTidy::assign_chains( const clipper::MiniMol& mol, const std::vector<int> seqnums, const double rmsd, const int nmin )
{
  const int nchn = mol.size();

  // source chain clash data
  const std::vector<int> num_seq = sequence_count( mol );
  const clipper::Array2d<int> seqflg = sequence_flags( mol );
  const int lseq = seqflg.cols();

  // build matrix of which chains superpose
  clipper::Array2d<int> super( nchn, nchn, false );
  for ( int c1 = 0; c1 < nchn-1; c1++ )
    for ( int c2 = c1 + 1; c2 < nchn; c2++ )
      if ( seqnums[c1] >= 0 && seqnums[c2] >= 0 &&
	   seqnums[c1] == seqnums[c2] ) {
	const clipper::RTop_orth rtop =
	  ProteinTools::superpose( mol[c2], mol[c1], rmsd, nmin, nmin );
	if ( ! rtop.is_null() ) super(c1,c2) = super(c2,c1) = true;
      }

  /*
  for ( int c1 = 0; c1 < nchn; c1++ ) std::cout << seqnums[c1] << ",";
  std::cout << std::endl;
  for ( int c1 = 0; c1 < nchn; c1++ ) {
    for ( int c2 = 0; c2 < nchn; c2++ ) std::cout << super(c1,c2);
    std::cout << std::endl;
  }
  */

  // make initial chain numbers list
  std::vector<int> chnnums( nchn, -1 );

  // find out how many sequences are present
  int nseqs = 0;
  for ( int c = 0; c < seqnums.size(); c++ )
    nseqs = std::max( nseqs, seqnums[c]+1 );

  // find the largest closed groups of ncs related chains for each sequence
  int nsrc = 0;
  for ( int s = 0; s < nseqs; s++ ) {
    std::vector<int> active( nchn ), used_best( nchn, false ), used;
    // make a list of chains matching this sequence 
    for ( int c = 0; c < nchn; c++ )
      active[c] = ( seqnums[c] == s );
    // find the biggest group of superposable chains
    best_closed_ncs_group( super, num_seq, active, used_best, used );
    // label them as starting groups
    for ( int c = 0; c < nchn; c++ )
      if ( used_best[c] ) chnnums[c] = nsrc++;
  }
  // std::cout << "SEQ: "; for ( int c = 0; c < nchn; c++ ) std::cout << clipper::String( seqnums[c], 5 ) << ","; std::cout << std::endl;
  // std::cout << "CNN: "; for ( int c = 0; c < nchn; c++ ) std::cout << clipper::String( chnnums[c], 5 ) << ","; std::cout << std::endl;

  // source chain map
  const double rad = 5.0;
  clipper::MAtomNonBond nb( mol, rad );

  // now expand chain numbers list
  for ( int cyc = 0; cyc < nchn; cyc++ ) {
    std::vector<double> scr1(nchn,0.0), scr2(nchn,0.0), scr(nchn,0.0);
    std::vector<int> chn(nchn,-1);

    // search over chains to find the best one to merge
    for ( int c = 0; c < nchn; c++ ) 
      if ( chnnums[c] < 0 && seqnums[c] >= 0 && num_seq[c] > 0 ) {
	std::vector<double> sscr1(nsrc,0.0), sscr2(nsrc,0.0), sscr(nsrc,0.0);

	// score by intimacy to known chains
	std::vector<int> refcount =
	  count_contacts( nb, mol, chnnums, mol[c], rad );
	for ( int s = 0; s < nsrc; s++ )
	  sscr1[s] = double(refcount[s]) / double(num_seq[c]);

	// score by overlap with known chains
	std::vector<int> clscount( nsrc, 0 );
	for ( int c1 = 0; c1 < nchn; c1++ )
	  if ( chnnums[c1] >= 0 )
	    for ( int r = 0; r < lseq; r++ )
	      if ( seqflg(c,r) && seqflg(c1,r) ) clscount[chnnums[c1]]++;
	for ( int s = 0; s < nsrc; s++ )
	  sscr2[s] = double(clscount[s]) / double(num_seq[c]);

	// filter out sources with clashing chain ids
	std::vector<int> sameseq( nsrc, true );
	for ( int c1 = 0; c1 < nchn; c1++ )
	  if ( chnnums[c1] >= 0 && seqnums[c1] != seqnums[c] )
	    sameseq[chnnums[c1]] = false;

	// total score
	for ( int s = 0; s < nsrc; s++ )
	  if ( sameseq[s] )
	    sscr[s] = sscr1[s] - 2.0*sscr2[s];

	// get best score for this chain
	int schn = 0;
	for ( int s = 0; s < nsrc; s++ )
	  if ( sscr[s] > sscr[schn] ) schn = s;
	chn[c] = schn;
	scr[c] = sscr[schn];
	scr1[c] = sscr1[schn];
	scr2[c] = sscr2[schn];
      }

    //std::cout << "CNN: "; for ( int c = 0; c < nchn; c++ ) std::cout << clipper::String( chnnums[c], 5 ) << ","; std::cout << std::endl;
    //std::cout << "CHN: "; for ( int c = 0; c < nchn; c++ ) std::cout << clipper::String( chn[c], 5 ) << ","; std::cout << std::endl;
    //std::cout << "SC1: "; for ( int c = 0; c < nchn; c++ ) std::cout << clipper::String( scr1[c], 5 ) << ","; std::cout << std::endl;
    //std::cout << "SC2: "; for ( int c = 0; c < nchn; c++ ) std::cout << clipper::String( scr2[c], 5 ) << ","; std::cout << std::endl << std::endl;
    // find best chain to assign
    int c = 0;
    for ( int c1 = 0; c1 < nchn; c1++ ) if ( scr[c1] > scr[c] ) c = c1;

    if ( scr[c] <= 0.0 ) break;

    // add new chain to source
    chnnums[c] = chn[c];
  }

  return chnnums;
}


std::vector<int> ModelTidy::chain_renumber( clipper::MiniMol& mol, const clipper::MMoleculeSequence& seq )
{
  std::vector<int> result( mol.size(), -1 );
  for ( int c = 0; c < mol.size(); c++ )
    result[c] = chain_renumber( mol[c], seq );
  return result;
}


int ModelTidy::chain_renumber( clipper::MPolymer& mp, const clipper::MMoleculeSequence& seq )
{
  // initial numbering
  for ( int i = 0; i < mp.size(); i++ ) mp[i].set_seqnum( i+1 );
  if ( seq.size() == 0 ) return -1;

  // convert sequences to unique strings
  clipper::String chnseq = ProteinTools::chain_sequence( mp );
  std::vector<clipper::String> seqs( seq.size() );
  for ( int chn = 0; chn < seq.size(); chn++ ) {
    clipper::String s = "";
    for ( int res = 0; res < seq[chn].sequence().length(); res++ )
      s += ProteinTools::residue_code_1(ProteinTools::residue_index_1(seq[chn].sequence().substr(res,1)));
    seqs[chn] = s;
  }
  
  // now find best match
  int bestchn = -1;
  int bestscr = 0;
  clipper::MSequenceAlign align( clipper::MSequenceAlign::LOCAL,
  				 1.0, -1000.0, -4.0 );
  std::pair<std::vector<int>,std::vector<int> > result;
  for ( int chn = 0; chn < seqs.size(); chn++ ) {
    const clipper::String& seqseq = seqs[chn];
    result = align( chnseq, seqseq );
    int scr = 0;
    for ( int i = 0; i < result.first.size(); i++ ) {
      if ( result.first[i] >= 0 ) {
	if ( chnseq[i] == seqseq[result.first[i]] )
 	  scr++;
  	else
  	  scr--;
      }
    }
    if ( scr > bestscr ) {
      bestchn = chn;
      bestscr = scr;
    }
  }
  if ( bestchn < 0 ) return -1;
  
  // now number residues
  clipper::String truseq = seqs[bestchn];
  result = align( chnseq, truseq );
  std::vector<int> nums = result.first;
  
  /*
  std::cout << bestchn << " " << bestscr << std::endl;
  for ( int i = 0; i < chnseq.size(); i++ )
    std::cout << chnseq[i];
  std::cout << std::endl;
  for ( int i = 0; i < chnseq.size(); i++ )
    std::cout << (nums[i]>=0&&nums[i]<truseq.length() ? truseq[nums[i]] : '-');
  std::cout << std::endl;
  */
  
  // find bounds of sequenced region
  int i0, i1, i;
  for ( i0 = 0;    i0 < nums.size(); i0++ ) if ( nums[i0] >= 0 ) break;
  for ( i1 = nums.size()-1; i1 >= 0; i1-- ) if ( nums[i1] >= 0 ) break;
  if ( i0 < nums.size() )
    for ( i = i0 - 1; i >= 0;          i-- ) nums[i] = nums[i+1] - 1;
  if ( i1 >= 0 )
    for ( i = i1 + 1; i < nums.size(); i++ ) nums[i] = nums[i-1] + 1;
  
  // renumber the model, with inscodes if necessary
  const clipper::String inscodes = " ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  const int l = inscodes.length();
  for ( i = 0; i < nums.size(); i++ ) nums[i] = nums[i] * l;
  for ( i = 1; i < nums.size(); i++ )
    if ( nums[i] <= nums[i-1] ) nums[i] = nums[i-1]+1;
  for ( i = 0; i < mp.size(); i++ ) {
    const int nres = (nums[i]/l) + 1;
    const int nins = clipper::Util::mod( nums[i], l );
    if ( nins == 0 ) mp[i].set_seqnum( nres );
    else             mp[i].set_seqnum( nres, inscodes.substr( nins, 1 ) );
  }

  // fill in sequence if reqiured
  for ( i = 0; i < mp.size(); i++ )
    if ( mp[i].type() == "###" ) {
      if ( nums[i] % l == 0 ) {
	const int s = nums[i]/l;
	const int t = ProteinTools::residue_index( truseq[s] );
	mp[i].set_type( ProteinTools::residue_code_3( t ) );
      } else {
	mp[i].set_type( "UNK" );
      }
    }

  /*
  for ( int i = 0; i < chnseq.size()-1; i++ )
    if ( nums[i+1] != nums[i]+1 ) std::cout << "! " << mp.id() << " " << nums[i] << " " << nums[i+1] << std::endl;
  */
  return bestchn;
}


bool ModelTidy::update_model( clipper::MiniMol& mol, const clipper::MiniMol& mol_mr, std::vector<int> chnnums ) const
{
  const clipper::Spacegroup& spgr = mol.spacegroup();
  const clipper::Cell&       cell = mol.cell();

  // current number of chains
  const int nchn = mol.size();
  // ideal number of chains (chnnums)
  int nsrc = 0;
  for ( int i = 0; i < chnnums.size(); i++ )
    nsrc = std::max( nsrc, chnnums[i]+1 );

  // chain ids
  std::vector<clipper::String> chnids;
  for ( int c = 0; c < mol_mr.size(); c++ ) chnids.push_back( mol_mr[c].id() );
  for ( int i = chnids.size(); i < nsrc; i++ ) chnids.push_back( "" );

  // assemble new chains
  clipper::MiniMol molnew( spgr, cell );
  for ( int c2 = 0; c2 < nsrc; c2++ ) {
    clipper::MPolymer mp;
    mp.set_id( chnids[c2] );

    // sort fragments by sequence number
    std::vector<std::pair<int,int> > srcsort;
    for ( int c1 = 0; c1 < nchn; c1++ ) {
      if ( chnnums[c1] == c2 ) {
	srcsort.push_back(std::pair<int,int>(mol[c1][0].seqnum(),c1));
      }
    }
    std::sort( srcsort.begin(), srcsort.end() );
    for ( int i = 0; i < srcsort.size(); i++ ) {
      int c1 = srcsort[i].second;
      for ( int r1 = 0; r1 < mol[c1].size(); r1++ )
	mp.insert( mol[c1][r1] );
      if ( verbose_ ) std::cout << "Adding chain " << mol[c1].id() << " to chain " << c2 << " " << mp.id() << std::endl;
    }

    // record chain gaps
    std::vector<bool> gaps( mp.size(), false );
    for ( int i = 1; i < gaps.size(); i++ )
      gaps[i] = ! clipper::MMonomer::protein_peptide_bond( mp[i-1], mp[i] );

    // renumber the model, with inscodes if necessary
    const clipper::String inscodes = " ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const int l = inscodes.length();
    std::vector<int> nums( mp.size() );
    for ( int i = 0; i < nums.size(); i++ ) nums[i] = mp[i].seqnum();
    for ( int i = 0; i < nums.size(); i++ ) nums[i] = nums[i] * l;
    for ( int i = 1; i < nums.size(); i++ )
      if ( gaps[i] ) nums[i] = std::max( nums[i], nums[i-1]-nums[i-1]%l+2*l );
      else           nums[i] = std::max( nums[i], nums[i-1]+1 );
    for ( int i = 0; i < mp.size(); i++ ) {
      const int nres = (nums[i]/l) + 1;
      const int nins = clipper::Util::mod( nums[i], l );
      if ( nins == 0 ) mp[i].set_seqnum( nres );
      else             mp[i].set_seqnum( nres, inscodes.substr( nins, 1 ) );
    }

    if ( mp.size() > 0 ) molnew.insert( mp );
  }

  // assemble unassigned chains
  clipper::MiniMol molunk( spgr, cell );
  for ( int c = 0; c < nchn; c++ )
    if ( chnnums[c] < 0 ) molunk.insert( mol[c] );
  if ( verbose_ ) std::cout << "Chains -  Old: " << nchn << "  Assigned: " << molnew.size() << "  Unassigned: " << molunk.size() << std::endl;
  // copy chains across
  for ( int c = 0; c < molunk.size(); c++ ) molunk[c].set_id( "" );
  for ( int c = 0; c < molunk.size(); c++ ) molnew.insert( molunk[c] );

  mol = molnew;
  return true;
}


bool ModelTidy::move_chain( clipper::MPolymer& mp1, const clipper::MPolymer& mp0, const clipper::Spacegroup& spgr, const clipper::Cell& cell )
{
  clipper::Coord_orth co;
  double no;

  // get COM of target chain
  co = clipper::Coord_orth( 0.0, 0.0, 0.0 );
  no = 0.0;
  for ( int r1 = 0; r1 < mp0.size(); r1++ )
    for ( int a1 = 0; a1 < mp0[r1].size(); a1++ ) {
      co += mp0[r1][a1].coord_orth();
      no += 1.0;
    }
  clipper::Coord_orth com0 = (1.0/no) * co;

  // get COM of moving chain
  co = clipper::Coord_orth( 0.0, 0.0, 0.0 );
  no = 0.0;
  for ( int r1 = 0; r1 < mp1.size(); r1++ )
    for ( int a1 = 0; a1 < mp1[r1].size(); a1++ ) {
      co += mp1[r1][a1].coord_orth();
      no += 1.0;
    }
  clipper::Coord_orth com1 = (1.0/no) * co;

  // find symop and lattice shift which brings it close to com0
  clipper::Coord_frac cf0 = com0.coord_frac( cell );
  clipper::Coord_frac cf1 = com1.coord_frac( cell );
  double r2min = 1.0e9;
  int symin = 0;
  clipper::Coord_frac cfmin( 0.0, 0.0, 0.0 );
  for ( int s = 0; s < spgr.num_symops(); s++ ) {
    clipper::Coord_frac cf2 = spgr.symop(s) * cf1;
    clipper::Coord_frac cf3 = cf2.lattice_copy_near( cf0 );
    double r2 = (cf3-cf0).lengthsq(cell);
    if ( r2 < r2min ) {
      r2min = r2;
      symin = s;
      cfmin = cf3 - cf2;
    }
  }

  // apply the shifts
  for ( int r1 = 0; r1 < mp1.size(); r1++ )
    for ( int a1 = 0; a1 < mp1[r1].size(); a1++ ) {
      cf1 = mp1[r1][a1].coord_orth().coord_frac( cell );
      cf1 = spgr.symop(symin) * cf1 + cfmin;
      mp1[r1][a1].set_coord_orth( cf1.coord_orth( cell ) );
    }

  return true;
}


std::vector<int> ModelTidy::count_contacts( const clipper::MAtomNonBond& nb, const clipper::MiniMol& mol, const std::vector<int>& chnnums, const clipper::MPolymer mp, const double& rad )
{
  // get number of sources
  int nsrc = 0;
  for ( int i = 0; i < chnnums.size(); i++ )
    nsrc = std::max( nsrc, chnnums[i]+1 );
  // now count contacts to sources
  std::vector<int> result( nsrc, 0 );
  // loop over target chain looking for residues
  for ( int r1 = 0; r1 < mp.size(); r1++ ) {
    if ( ProteinTools::residue_index_3( mp[r1].type() ) >= 0 ) {
      int a1 = mp[r1].lookup( " CA ", clipper::MM::ANY );
      if ( a1 >= 0 ) {
	// find nearby source model atoms
        const std::vector<clipper::MAtomIndexSymmetry> atoms =
          nb( mp[r1][a1].coord_orth(), rad );
        int chnn1 = -1;
	// loop over source model atoms
        for ( int i = 0; i < atoms.size(); i++ ) {
          const int c2 = atoms[i].polymer();
          const int r2 = atoms[i].monomer();
	  int chnn2 = chnnums[c2];
	  // if this atom has an index number then consider it
	  if ( chnn2 >= 0 ) {
	    if ( ProteinTools::residue_index_3( mol[c2][r2].type() ) >= 0 ) {
	      if ( chnn1 < 0 ) {
		chnn1 = chnn2;
	      } else if ( chnn1 != chnn2 ) {
		chnn1 = -1;
		break;
	      }
	    }
	  }
        }
        if ( chnn1 >= 0 ) result[chnn1]++;
      }
    }
  }
  return result;
}


std::vector<int> ModelTidy::sequence_count( const clipper::MiniMol& mol )
{
  std::vector<int> result( mol.size(), 0 );
  for ( int c = 0; c < mol.size(); c++ )
    for ( int r = 0; r < mol[c].size(); r++ )
      if ( ProteinTools::residue_index_3( mol[c][r].type() ) >= 0 )
	result[c]++;
  return result;
}


clipper::Array2d<int> ModelTidy::sequence_flags( const clipper::MiniMol& mol )
{
  const int nchn = mol.size();

  // make list of used residue numbers
  int smin(1000000), smax(-1000000);
  for ( int c = 0; c < mol.size(); c++ )
    for ( int r = 0; r < mol[c].size(); r++ )
      if ( ProteinTools::residue_index_3( mol[c][r].type() ) >= 0 ) {
	int s = mol[c][r].seqnum();
	smin = std::min( smin, s );
	smax = std::max( smax, s );
      }
  int nseq = smax - smin + 1;
  if ( nseq <= 0 ) return clipper::Array2d<int>( 0, 0 );

  // make flags of used residue numbers
  clipper::Array2d<int> seqflg( nchn, nseq, false );
  for ( int c = 0; c < mol.size(); c++ )
    for ( int r = 0; r < mol[c].size(); r++ )
      if ( ProteinTools::residue_index_3( mol[c][r].type() ) >= 0 ) {
	int s = mol[c][r].seqnum();
	seqflg(c,s-smin) = true;
      }

  return seqflg;
}


void ModelTidy::best_closed_ncs_group( const clipper::Array2d<int>& super, const std::vector<int>& num_seq, const std::vector<int>& active, std::vector<int>& used_best, std::vector<int>& used )
{
  const int nchn = used_best.size();
  const int chn  = used.size();
  if ( chn == nchn ) {
    // finished recursion - test if better solution found
    int n1(0),n2(0),r1(0),r2(0);
    for ( int c = 0; c < nchn; c++ ) {
      if ( used[c] )      { n1++; r1+=num_seq[c]; }
      if ( used_best[c] ) { n2++; r2+=num_seq[c]; }
    }
    if ( n1 > n2 || ( n1 == n2 && r1 > r2 ) ) used_best = used;
    //for ( int c = 0; c < nchn; c++ ) std::cout << used[c];
    //std::cout << " " << n1 << " " << r1 << std::endl;
  } else {
    // check if we can use this chain
    bool can_use = active[chn];
    for ( int c = 0; c < chn; c++ )
      if ( used[c] ) can_use = can_use && super(c,chn);
    if ( can_use ) {
      // if so, try using it
      used.push_back( true );
      best_closed_ncs_group( super, num_seq, active, used_best, used );
      used.pop_back();
    }
    // try not using it
    used.push_back( false );
    best_closed_ncs_group( super, num_seq, active, used_best, used );
    used.pop_back();
  }
}
