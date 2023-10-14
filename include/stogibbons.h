#ifndef STOGIBBONS_H
#define STOGIBBONS_H

#include "estimator.h"
#include <set>

/**
 * Used below in BarYossef function
 */
template <class T1, class T2> class HashValuePair {
public:
  HashValuePair(T1 t1, const T2 &t2) : mhash(t1), mvalue(t2) {}
  bool operator==(const HashValuePair<T1, T2> &y) const {
    return mvalue == y.mvalue;
  }
  bool operator<(const HashValuePair<T1, T2> &y) const {
    return mhash < y.mhash;
  }
  T1 hash() const { return mhash; }
  T2 &value() const { return mvalue; }

private:
  T1 mhash;
  T2 mvalue;
};

/**
 * Here I decided to implement an algorithm by Bar-Yossef et al.
 * from the paper "counting distinct elements in a data stream".
 * On paper, it looks worse that Gibbons, but who knows how well it may fare?
 * Its implementation is remarkably easy.
 *  See http://www.ee.technion.ac.il/people/zivby/papers/f0/f0.ps
 *
 */
template <class T>
uint64 BarYossef(istream &in, vector<string> &groupby, int &N, int L = 19,
                 uint storagelimit = 512, char delimiter = ',',
                 int rng = rand_rng) {
  // int nbofbitmap = 32;
  set<HashValuePair<uint64, vector<T>>> buffer;

  const uint n = groupby.size(); // number of dimensions
  uint64 largesthash(0);
  StringConverter<T> sc;

  TupleFactHash<T> cph(n, L, rng);
  string line; // Read a line of data
  vector<T> attribute_value(n);
  vector<string> attribute_list; // Pour stocker les attributs de dimension et
                                 // leur valeur de chaque faits
  vector<unsigned int> dimindice; // Pour stocker les indices des colonnes
                                  // correspondant aux attributs du group by

  N = 0;
  // assert(buffer.size() == 0);  // if not, we needed to clear...

  // Récupérer les indices des colones des dimensions
  getline(in, line); // La première ligne du fichier de données contient les
                     // noms des dimensions

  string delimitertostring(1, delimiter);
  split(line, attribute_list, delimitertostring);

  for (unsigned int j = 0; j < groupby.size(); ++j)
    for (unsigned int i = 0; i < attribute_list.size(); ++i)
      if (groupby[j] == attribute_list[i]) {
        dimindice.push_back(i);
        break;
      }
  sort(dimindice.begin(),
       dimindice.end()); // Ascending sort of the dimension indices

  assert(dimindice.size() == n);
  char chline[256];
  while (in.getline(chline, MAXLENGTHOFALINE, delimiter)) {
    ++N;

    uint col(0);
    bool endline(false);
    for (uint k = 0; k < dimindice.size(); ++k) {
      while (col < dimindice[k]) {
        ++col;
        if (col == attribute_list.size() - 1) {
          in.getline(chline, MAXLENGTHOFALINE);
          endline = true;
        } else
          in.getline(chline, MAXLENGTHOFALINE, delimiter);
      }
      attribute_value[k] = sc.toType(chline);
    }
    if (!endline)
      in.getline(chline, MAXLENGTHOFALINE);
    uint64 myhash = cph.slow(attribute_value); // hashing
    if (buffer.size() < storagelimit) {        // buffer not yet initialized
      buffer.insert(HashValuePair<uint64, vector<T>>(myhash, attribute_value));
      largesthash = largesthash > myhash ? largesthash : myhash;
    } else if (myhash < largesthash) { // we only add if it is worth it
      buffer.insert(HashValuePair<uint64, vector<T>>(myhash, attribute_value));
      if (buffer.size() > storagelimit) {
        buffer.erase(--buffer.end());
        largesthash = (--buffer.end())->hash();
      }
    }
  }
  uint64 lastvalue = (--buffer.end())->hash();
  return (uint64)round(buffer.size() * pow(2.0, L) / lastvalue);
}

struct adder  {
  adder() : sum(0) {}
  uint64 sum;
  void operator()(uint64 x) { sum += x; }
};


template <class T>
uint64 StoGibbonsTirthapura(istream &in, vector<string> &groupby, int &N,
                            int L = 19, int storagelimit = 512,
                            int nbofbitmap = 32, char delimiter = ',',
                            int rng = rand_rng) {
  // int nbofbitmap = 32;
  vector<map<vector<T>, uint64>> buffer(nbofbitmap);

  const uint n = groupby.size(); // number of dimensions
  StringConverter<T> sc;

  // Estimating simple functions on the union of data streams by Gibbons and
  // Tirthapura
  TupleFactHash<T> cph(n, L, rng);
  vector<uint> t(nbofbitmap, 0);
  vector<uint> mask(nbofbitmap,
                    0); // = (1 << (t)) - 1;   // which is initially 0
  string line;          // Read a line of data
  vector<T> attribute_value(n);
  vector<string> attribute_list; // Pour stocker les attributs de dimension et
                                 // leur valeur de chaque faits
  vector<unsigned int> dimindice; // Pour stocker les indices des colonnes
                                  // correspondant aux attributs du group by

  N = 0;
  // assert(buffer.size() == 0);  // if not, we needed to clear...

  // Récupérer les indices des colones des dimensions
  getline(in, line); // La première ligne du fichier de données contient les
                     // noms des dimensions

  string delimitertostring(1, delimiter);
  split(line, attribute_list, delimitertostring);

  for (unsigned int j = 0; j < groupby.size(); ++j)
    for (unsigned int i = 0; i < attribute_list.size(); ++i)
      if (groupby[j] == attribute_list[i]) {
        dimindice.push_back(i);
        break;
      }
  sort(dimindice.begin(),
       dimindice.end()); // Ascending sort of the dimension indices

  assert(dimindice.size() == n);
  char chline[256];
  while (in.getline(chline, MAXLENGTHOFALINE, delimiter)) {
    ++N;

    uint col(0);
    bool endline(false);
    for (uint k = 0; k < dimindice.size(); ++k) {
      while (col < dimindice[k]) {
        ++col;
        if (col == attribute_list.size() - 1) {
          in.getline(chline, MAXLENGTHOFALINE);
          endline = true;
        } else
          in.getline(chline, MAXLENGTHOFALINE, delimiter);
      }
      attribute_value[k] = sc.toType(chline);
    }
    if (!endline)
      in.getline(chline, MAXLENGTHOFALINE);

    uint myhash = cph.slow(attribute_value); // hashing
    int alpha = myhash % nbofbitmap;
    uint oldhash = myhash / nbofbitmap;

    // buffer[attribute_value] = (1LL << 32) | trailingzeros(oldhash,L);  //
    // count= 1
    uint projection = oldhash & mask[alpha];
    if (projection == 0) {
      typename map<vector<T>, uint64>::iterator it =
          buffer[alpha].find(attribute_value);
      if (it == buffer[alpha].end()) {
        // buffer[attribute_value] = (1LL << 32) | trailingzeros(oldhash,L);
        buffer[alpha][attribute_value] =
            trailingzeros(oldhash, L); //(1LL << 32) |
        while ((buffer[alpha].size() > (uint)storagelimit / nbofbitmap)) {
          t[alpha] += 1;
          mask[alpha] = (1 << (t[alpha])) - 1;
          assert(t[alpha] < (uint)L);
          for (typename map<vector<T>, uint64>::iterator j =
                   buffer[alpha].begin();
               j != buffer[alpha].end();)
            // if( (j->second & 0xffffffff)  < t) buffer.erase(j++);
            if ((j->second) < t[alpha])
              buffer[alpha].erase(j++); // & 0xffffffff
            else
              ++j;
        }
      } // else // previously seen, so bump count by 1
      // buffer[attribute_value] += (1LL << 32);
    }
  }
  vector<int> estimates(nbofbitmap, 0);
  for (int k = 0; k < nbofbitmap; ++k)
    estimates[k] = buffer[k].size() * (1 << t[k]) * nbofbitmap;
  sort(estimates.begin(), estimates.end());
  // int estimate = 0;
  if (nbofbitmap < 1)
    cerr << "oops! expect garbage!!!" << endl;
  if (nbofbitmap == 1)
    return estimates[0];
  if (nbofbitmap % 2 == 0)
    return (int)round(
        (estimates[nbofbitmap / 2 - 1] + estimates[(nbofbitmap / 2)]) / 2.0);
  return estimates[nbofbitmap / 2];
}
template <class T>
uint64 StoSumGibbonsTirthapura(istream &in, vector<string> &groupby, int &N,
                               int L = 19, int storagelimit = 512,
                               int nbofbitmap = 32, char delimiter = ',',
                               int rng = rand_rng) {
  // int nbofbitmap = 32;
  vector<map<vector<T>, uint64>> buffer(nbofbitmap);

  const uint n = groupby.size(); // number of dimensions
  StringConverter<T> sc;

  // Estimating simple functions on the union of data streams by Gibbons and
  // Tirthapura
  TupleFactHash<T> cph(n, L, rng);
  vector<uint> t(nbofbitmap, 0);
  vector<uint> mask(nbofbitmap,
                    0); // = (1 << (t)) - 1;   // which is initially 0
  string line;          // Read a line of data
  vector<T> attribute_value(n);
  vector<string> attribute_list; // Pour stocker les attributs de dimension et
                                 // leur valeur de chaque faits
  vector<unsigned int> dimindice; // Pour stocker les indices des colonnes
                                  // correspondant aux attributs du group by

  N = 0;
  // assert(buffer.size() == 0);  // if not, we needed to clear...

  // Récupérer les indices des colones des dimensions
  getline(in, line); // La première ligne du fichier de données contient les
                     // noms des dimensions

  string delimitertostring(1, delimiter);
  split(line, attribute_list, delimitertostring);

  for (unsigned int j = 0; j < groupby.size(); ++j)
    for (unsigned int i = 0; i < attribute_list.size(); ++i)
      if (groupby[j] == attribute_list[i]) {
        dimindice.push_back(i);
        break;
      }
  sort(dimindice.begin(),
       dimindice.end()); // Ascending sort of the dimension indices

  assert(dimindice.size() == n);
  char chline[256];
  while (in.getline(chline, MAXLENGTHOFALINE, delimiter)) {
    ++N;

    uint col(0);
    bool endline(false);
    for (uint k = 0; k < dimindice.size(); ++k) {
      while (col < dimindice[k]) {
        ++col;
        if (col == attribute_list.size() - 1) {
          in.getline(chline, MAXLENGTHOFALINE);
          endline = true;
        } else
          in.getline(chline, MAXLENGTHOFALINE, delimiter);
      }
      attribute_value[k] = sc.toType(chline);
    }
    if (!endline)
      in.getline(chline, MAXLENGTHOFALINE);

    uint myhash = cph.slow(attribute_value); // hashing
    int alpha = myhash % nbofbitmap;
    uint oldhash = myhash / nbofbitmap;

    // buffer[attribute_value] = (1LL << 32) | trailingzeros(oldhash,L);  //
    // count= 1
    uint projection = oldhash & mask[alpha];
    if (projection == 0) {
      typename map<vector<T>, uint64>::iterator it =
          buffer[alpha].find(attribute_value);
      if (it == buffer[alpha].end()) {
        // buffer[attribute_value] = (1LL << 32) | trailingzeros(oldhash,L);
        buffer[alpha][attribute_value] =
            trailingzeros(oldhash, L); //(1LL << 32) |
        while ((buffer[alpha].size() > (uint)storagelimit / nbofbitmap)) {
          t[alpha] += 1;
          mask[alpha] = (1 << (t[alpha])) - 1;
          assert(t[alpha] < (uint)L);
          for (typename map<vector<T>, uint64>::iterator j =
                   buffer[alpha].begin();
               j != buffer[alpha].end();)
            // if( (j->second & 0xffffffff)  < t) buffer.erase(j++);
            if ((j->second) < t[alpha])
              buffer[alpha].erase(j++); // & 0xffffffff
            else
              ++j;
        }
      } // else // previously seen, so bump count by 1
      // buffer[attribute_value] += (1LL << 32);
    }
  }
  vector<int> estimates(nbofbitmap, 0);
  for (int k = 0; k < nbofbitmap; ++k) {
    estimates[k] = buffer[k].size() * (1 << t[k]);
    if (stoverbose)
      cout << buffer[k].size() << " " << t[k] << endl;
  }
  return for_each(estimates.begin(), estimates.end(), adder()).sum;
}

/**
 * possibly garbage.
 */
template <class T>
uint64
StoWeightedSumGibbonsTirthapura(istream &in, vector<string> &groupby, int &N,
                                int L = 19, int storagelimit = 512,
                                int nbofbitmap = 32, char delimiter = ',',
                                int rng = rand_rng) {
  // int nbofbitmap = 32;
  vector<map<vector<T>, uint64>> buffer(nbofbitmap);

  const uint n = groupby.size(); // number of dimensions
  StringConverter<T> sc;

  // Estimating simple functions on the union of data streams by Gibbons and
  // Tirthapura
  TupleFactHash<T> cph(n, L, rng);
  vector<uint> t(nbofbitmap, 0);
  vector<uint> mask(nbofbitmap,
                    0); // = (1 << (t)) - 1;   // which is initially 0
  string line;          // Read a line of data
  vector<T> attribute_value(n);
  vector<string> attribute_list; // Pour stocker les attributs de dimension et
                                 // leur valeur de chaque faits
  vector<unsigned int> dimindice; // Pour stocker les indices des colonnes
                                  // correspondant aux attributs du group by

  N = 0;
  // assert(buffer.size() == 0);  // if not, we needed to clear...

  // Récupérer les indices des colones des dimensions
  getline(in, line); // La première ligne du fichier de données contient les
                     // noms des dimensions

  string delimitertostring(1, delimiter);
  split(line, attribute_list, delimitertostring);

  for (unsigned int j = 0; j < groupby.size(); ++j)
    for (unsigned int i = 0; i < attribute_list.size(); ++i)
      if (groupby[j] == attribute_list[i]) {
        dimindice.push_back(i);
        break;
      }
  sort(dimindice.begin(),
       dimindice.end()); // Ascending sort of the dimension indices

  assert(dimindice.size() == n);
  char chline[256];
  while (in.getline(chline, MAXLENGTHOFALINE, delimiter)) {
    ++N;

    uint col(0);
    bool endline(false);
    for (uint k = 0; k < dimindice.size(); ++k) {
      while (col < dimindice[k]) {
        ++col;
        if (col == attribute_list.size() - 1) {
          in.getline(chline, MAXLENGTHOFALINE);
          endline = true;
        } else
          in.getline(chline, MAXLENGTHOFALINE, delimiter);
      }
      attribute_value[k] = sc.toType(chline);
    }
    if (!endline)
      in.getline(chline, MAXLENGTHOFALINE);

    uint myhash = cph.slow(attribute_value); // hashing
    int alpha = myhash % nbofbitmap;
    uint oldhash = myhash / nbofbitmap;

    uint projection = oldhash & mask[alpha];

    if (projection == 0) {
      typename map<vector<T>, uint64>::iterator it =
          buffer[alpha].find(attribute_value);
      if (it == buffer[alpha].end()) {
        // buffer[attribute_value] = (1LL << 32) | trailingzeros(oldhash,L);
        buffer[alpha][attribute_value] =
            trailingzeros(oldhash, L); //(1LL << 32) |
        while ((buffer[alpha].size() > (uint)storagelimit / nbofbitmap)) {
          t[alpha] += 1;
          mask[alpha] = (1 << (t[alpha])) - 1;
          assert(t[alpha] < (uint)L);
          for (typename map<vector<T>, uint64>::iterator j =
                   buffer[alpha].begin();
               j != buffer[alpha].end();)
            // if( (j->second & 0xffffffff)  < t) buffer.erase(j++);
            if ((j->second) < t[alpha])
              buffer[alpha].erase(j++); // & 0xffffffff
            else
              ++j;
        }
      } // else // previously seen, so bump count by 1
      // buffer[attribute_value] += (1LL << 32);
    }
  }
  vector<int> counts(nbofbitmap, 0);
  // vector<int> counts(nbofbitmap,0);
  for (int k = 0; k < nbofbitmap; ++k) {
    counts[k] = buffer[k].size();
    // unrolling, this is a heuristic daniel thought of
    if (counts[k] == 0) {
      t[k] -= 1;
      counts[k] += ((uint)storagelimit / nbofbitmap) + 1;
    }
  }
  double totalcount = for_each(counts.begin(), counts.end(), adder()).sum;
  double totalt = for_each(t.begin(), t.end(), adder()).sum;
  float phi = 1.0;
  if (storagelimit / nbofbitmap == 0) // Flajolet's case
    phi = 0.77351;
  if (stoverbose)
    cout << "nbofbitmap" << nbofbitmap << endl;
  if (stoverbose)
    cout << "totalt" << totalt << endl;
  return (uint)((totalcount / phi) * ((pow(2, (totalt / nbofbitmap)))));
}

#endif
