
//  $Id: front_end.H,v 1.5 2009/09/14 21:17:41 stanchen Exp $

/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
 *   @file front_end.H
 *   @brief Contains FrontEnd class for signal processing.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _FRONT_END_H
#define _FRONT_END_H

#include "util.H"

/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
 *   Front end class; does signal processing.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
class FrontEnd {
 public:
  /** Ctor; get configuration parameters from @p params. **/
  FrontEnd(const map<string, string>& params = ParamsType())
      : params_(params) {}

  /** Process input signal from @p inAudio; place acoustic
   *   feature vectors in @p out_feats.
   **/
  void get_feats(const matrix<double>& inAudio, matrix<double>& out_feats) const;

 private:
  /** Do windowing. **/
  void do_window(const matrix<double>& in_feats, matrix<double>& out_feats) const;

  /** Do FFT. **/
  void do_fft(const matrix<double>& in_feats, matrix<double>& out_feats) const;

  /** Do Mel binning. **/
  void do_melbin(const matrix<double>& in_feats, matrix<double>& out_feats) const;

  /** Do discrete cosine transform. **/
  void do_dct(const matrix<double>& in_feats, matrix<double>& out_feats) const;

 private:
  /** Stores copy of input parameters. **/
  map<string, string> params_;
};


struct Mel_bin_height {
    Mel_bin_height() {
        start_dim_idx = -1;
    }
    bool is_empty() {
        return start_dim_idx == -1;
    }
    int start_dim_idx;
    vector <double> heights;
};

void add_bin_height(vector<Mel_bin_height>& mel_bin_heights,
                    int const bin_idx,
                    int const dim_idx,
                    double const height);


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#endif
