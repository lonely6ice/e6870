
//  $Id: front_end.C,v 1.2 2016/01/23 03:15:23 stanchen Exp $

#include "front_end.H"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/** Module for doing windowing. **/
void FrontEnd::do_window(const matrix<double>& in_feats,
                         matrix<double>& out_feats) const {
  //  Get parameters.
  //  Input samples per second.
  double sample_rate = get_float_param(params_, "window.sample_rate", 20000.0);
  //  Output frames per second.
  double frames_per_sec =
      get_float_param(params_, "window.frames_per_sec", 100.0);
  //  Width of each window, in seconds.
  double window_width = get_float_param(params_, "window.window_size", 0.025);
  //  Whether to do Hamming or rectangular windowing.
  bool do_Hamming = get_bool_param(params_, "window.hamming", true);

  //  Get number of input samples.
  int in_samp_cnt = in_feats.size1();
  if (in_feats.size2() != 1)
    throw runtime_error("Windowing expected vector input.");

  //  Input sampling period in seconds.
  double sample_period = 1.0 / sample_rate;
  //  Output frame period, in seconds.
  double frame_period = 1.0 / frames_per_sec;
  //  Number of samples per window.
  int samp_per_window = (int)(window_width / sample_period + 0.5);
  //  Number of samples to shift between each window.
  int samp_shift = (int)(frame_period / sample_period + 0.5);
  //  Number of output frames.
  int out_frame_cnt = (in_samp_cnt - samp_per_window) / samp_shift + 1;

  //  Allocate output matrix and fill with zeros.
  out_feats.resize(out_frame_cnt, samp_per_window);
  out_feats.clear();

  //  BEGIN_LAB
  //
  //  Input:
  //      "in_feats", a matrix containing a single column holding the
  //      input samples for an utterance.  Each row holds a single sample.
  //
  //      in_feats(0 .. (in_samp_cnt - 1), 0)
  //
  //  Output:
  //      "out_feats", which should contain the result of windowing.
  //
  //      out_feats(0 .. (out_frame_cnt - 1), 0 .. (samp_per_window - 1))
  //
  //      Each row corresponds to a frame, and should hold the
  //      windowed samples for that frame.
  //      It has already been allocated to be of the correct size.
  //      If the boolean "do_Hamming" is true, then a Hamming
  //      window should be applied; otherwise, a rectangular
  //      window should be used.
  //
  //  See "in_samp_cnt", "samp_per_window", "samp_shift", and "out_frame_cnt"
  //  above for quantities you may (or may not) need for this computation.
  //
  //  When accessing matrices such as "in_feats" and "out_feats",
  //  use a syntax like "in_feats(frm_idx, dim_idx)" to access elements;
  //  using square brackets as in normal C arrays won't work.

  // cout << format("samp_shift %d\n") % samp_shift;
  // cout << format("out_frame_cnt %d\n") % out_frame_cnt;
  // cout << format("samp_perWindw %d\n") % samp_per_window;
  // cout << format("in_samp_cnt %d\n") % in_samp_cnt;
  if (do_Hamming) {
    // Hamming windows
    vector<double> hamming_window(samp_per_window, 0.0);
    double const PI = acos(-1.0);
    for (int i = 0; i < samp_per_window; ++i) {
      hamming_window[i] = 0.54 - 0.46 * cos(2 * PI * i / (samp_per_window - 1));
    }
    int frm_idx = 0;
    for (int samp_beg = 0; frm_idx < out_frame_cnt; samp_beg += samp_shift) {
      for (int i = 0; i < samp_per_window; ++i) {
        out_feats(frm_idx, i) = in_feats(samp_beg + i, 0) * hamming_window[i];
      }
      ++frm_idx;
    }
  } else {
    // Rectangular window
    int frm_idx = -1;
    for (int samp_beg = 0; frm_idx < out_frame_cnt; samp_beg += samp_shift) {
      copy(&in_feats(samp_beg, 0), &in_feats(samp_beg, 0) + samp_per_window, \
        &out_feats(++frm_idx, 0));
    }
  }
  //  END_LAB
}

/** Module for doing FFT. **/
void FrontEnd::do_fft(const matrix<double>& in_feats,
                      matrix<double>& out_feats) const {
  //  Make output dimension the smallest power of 2 at least as
  //  large as input dimension.
  int in_frame_cnt = in_feats.size1();
  int in_dim_cnt = in_feats.size2();
  int out_dim_cnt = 2;
  while (out_dim_cnt < in_dim_cnt) out_dim_cnt *= 2;

  //  Allocate output matrix and fill with zeros.
  out_feats.resize(in_frame_cnt, out_dim_cnt);
  out_feats.clear();

  //  Input:
  //      "in_feats", a matrix with each row holding the windowed
  //      values for that frame.
  //
  //      in_feats(0 .. (in_frame_cnt - 1), 0 .. (in_dim_cnt - 1))
  //
  //  Output:
  //      "out_feats", where an FFT should be applied to each
  //      row/frame of "in_feats".
  //
  //      out_feats(0 .. (in_frame_cnt - 1), 0 .. (out_dim_cnt - 1))
  //
  //      For a given row/frame "frm_idx", the real and imaginary
  //      parts of the FFT value for frequency i/(out_dim_cnt*T)
  //      where T is the sample period are held in
  //      out_feats(frm_idx, 2*i) and out_feats(frm_idx, 2*i+1),
  //      respectively.

  vector<double> fft_buf;
  for (int frm_idx = 0; frm_idx < in_frame_cnt; ++frm_idx) {
    copy_matrix_row_to_vector(in_feats, frm_idx, fft_buf);
    //  Pad window with zeros, if needed.
    fft_buf.resize(out_dim_cnt, 0.0);
    real_fft(fft_buf);
    copy_vector_to_matrix_row(fft_buf, out_feats, frm_idx);
  }
}

/** Module for mel binning. **/
// change to google name style
void FrontEnd::do_melbin(const matrix<double>& in_feats,
                         matrix<double>& out_feats) const {
  //  Number of mel bins to make.
  int num_bins = get_int_param(params_, "melbin.bins", 26);
  //  Whether to take log of output or not.
  bool do_log = get_bool_param(params_, "melbin.log", true);
  //  Input samples per second.
  double sample_rate = get_float_param(params_, "window.sample_rate", 20000.0);
  double sample_period = 1.0 / sample_rate;

  //  Retrieve number of frames and dimension of input feature vectors.
  int in_frame_cnt = in_feats.size1();
  int in_dim_cnt = in_feats.size2();
  int out_dim_cnt = num_bins;

  //  Allocate output matrix and fill with zeros.
  out_feats.resize(in_frame_cnt, out_dim_cnt);
  out_feats.clear();

  //  BEGIN_LAB
  //
  //  Input:
  //      "in_feats", holding the output of a real FFT.
  //
  //      in_feats(0 .. (in_frame_cnt - 1), 0 .. (in_dim_cnt - 1))
  //
  //  Output:
  //      "out_feats", which should contain the result of
  //      mel-binning.
  //
  //      out_feats(0 .. (in_frame_cnt - 1), 0 .. (out_dim_cnt - 1))
  //
  //      If the boolean "doLog" is true,
  //      then each value should be replaced with its natural
  //      logarithm, or 0 if its logarithm is negative.
  //      "out_feats" has been allocated to be of the correct size.
  //
  //  See "in_frame_cnt", "in_dim_cnt", "out_dim_cnt", and "sample_period"
  //  above for quantities you will need for this computation.

  // cout << format("sample_period %ls\n") % sample_period;
  // cout << format("input dim %d\n") % in_dim_cnt;

  // TODO: Optimize the compuatation order to avoid duplicate computation

  #define Mel(f) (1127 * log(1 + (f) / 700.0))
  // determine bin boundary
  vector<double> mel_bin_corners(out_dim_cnt + 2, 0.0);
  double const mel_f_max = Mel(1.0 / (2.0 * sample_period));
  for (int cor_idx = 0; cor_idx < out_dim_cnt + 2; ++cor_idx) {
      mel_bin_corners[cor_idx] = mel_f_max * cor_idx / (out_dim_cnt + 1);
  }
  double const mel_bin_with = mel_bin_corners[1];

  // determine bin hight
  vector<Mel_bin_height> mel_bin_heights(out_dim_cnt);
  int bin_idx = 0;
  for (int dim_idx = 0; dim_idx < in_dim_cnt / 2; ++dim_idx) {
    double mel_f = Mel(dim_idx / (in_dim_cnt * sample_period));
    if (mel_f > mel_bin_corners[bin_idx + 1]) {
      // find which half bin(in rising) include mel_f
      while (mel_f > mel_bin_corners[++bin_idx + 1]) {}
    }
    // cur bin weight in mel_f
    double height = (mel_f - mel_bin_corners[bin_idx]) / mel_bin_with;
    // add to current bin(in rising)
    add_bin_height(mel_bin_heights, bin_idx, dim_idx, height);
    // add to last bin(in falling)
    add_bin_height(mel_bin_heights, bin_idx - 1, dim_idx, 1 - height);
  }

  // calculate mel binning
  for (int frm_idx = 0; frm_idx < in_frame_cnt; ++frm_idx) {
    vector<double> amplitudes(in_dim_cnt / 2);
    for (int dim_idx = 0; dim_idx < in_dim_cnt / 2; ++dim_idx) {
      amplitudes[dim_idx] = sqrt(pow(in_feats(frm_idx, dim_idx * 2), 2) \
                               + pow(in_feats(frm_idx, dim_idx * 2 + 1), 2));
    }

    for (int bin_idx = 0; bin_idx < out_dim_cnt; ++bin_idx) {
      if (mel_bin_heights[bin_idx].is_empty()) {
        continue;
      }
      double sum = 0.0;
      for (int i = 0; i < int(mel_bin_heights[bin_idx].heights.size()); ++i) {
        sum += (amplitudes[mel_bin_heights[bin_idx].start_dim_idx + i] \
               * mel_bin_heights[bin_idx].heights[i]);
      }
      if (do_log) {
        sum = log(sum);
      }
      out_feats(frm_idx, bin_idx) = sum;
    }
  }
  //  END_LAB
}

/** Module for doing discrete cosine transform. **/
// change to google name style
void FrontEnd::do_dct(const matrix<double>& in_feats,
                      matrix<double>& out_feats) const {
  //  Number of DCT coefficients to output.
  int num_coeffs = get_int_param(params_, "dct.coeffs", 12);
  int in_frame_cnt = in_feats.size1();
  int in_dim_cnt = in_feats.size2();
  int out_dim_cnt = num_coeffs;

  //  Allocate output matrix and fill with zeros.
  out_feats.resize(in_frame_cnt, out_dim_cnt);
  out_feats.clear();

  //  BEGIN_LAB
  //
  //  Input:
  //      The matrix "in_feats", holding the output of mel-binning.
  //
  //      in_feats(0 .. (in_frame_cnt - 1), 0 .. (in_dim_cnt - 1))
  //
  //  Output:
  //      The matrix "out_feats", which should contain the result of
  //      applying the DCT.
  //
  //      out_feats(0 .. (in_frame_cnt - 1), 0 .. (out_dim_cnt - 1))
  //
  //      "out_feats" has been allocated to be of the correct size.
  //
  //  See "in_frame_cnt", "in_dim_cnt", and "out_dim_cnt" above
  //  for quantities you will need for this computation.

  double const PI = acos(-1.0);
  double const coeff = sqrt(2.0 / in_dim_cnt);
  for (int frm_idx = 0; frm_idx < in_frame_cnt; ++frm_idx) {
    for (int out_dim_idx = 0; out_dim_idx < out_dim_cnt; ++out_dim_idx) {
      double sum = 0.0;
      for (int in_dim_idx = 0; in_dim_idx < in_dim_cnt; ++in_dim_idx) {
        sum += in_feats(frm_idx, in_dim_idx) * \
               cos(PI * (out_dim_idx + 1) * (in_dim_idx + 0.5) / in_dim_cnt);
      }
      out_feats(frm_idx, out_dim_idx) = coeff * sum;
    }
  }

  //  END_LAB
}

/** Main signal processing routine.
 *   Calls each signal processing module in turn, unless
 *   parameter says not to.
 **/
void FrontEnd::get_feats(const matrix<double>& inAudio,
                         matrix<double>& out_feats) const {
  if (get_bool_param(params_, "frontend.null", false)) {
    out_feats = inAudio;
    return;
  }
  matrix<double> curFeats(inAudio);
  if (get_bool_param(params_, "frontend.window", true)) {
    do_window(curFeats, out_feats);
    out_feats.swap(curFeats);
  }
  if (get_bool_param(params_, "frontend.fft", true)) {
    do_fft(curFeats, out_feats);
    out_feats.swap(curFeats);
  }
  if (get_bool_param(params_, "frontend.melbin", true)) {
    do_melbin(curFeats, out_feats);
    out_feats.swap(curFeats);
  }
  if (get_bool_param(params_, "frontend.dct", true)) {
    do_dct(curFeats, out_feats);
    out_feats.swap(curFeats);
  }
  out_feats.swap(curFeats);
}

void add_bin_height(vector<Mel_bin_height>& mel_bin_heights,
                    int const bin_idx,
                    int const dim_idx,
                    double const height) {
    if (bin_idx >= 0 && bin_idx < int(mel_bin_heights.size())) {
        if (mel_bin_heights[bin_idx].is_empty()) {
            mel_bin_heights[bin_idx].start_dim_idx = dim_idx;
        }
        mel_bin_heights[bin_idx].heights.push_back(height);
    }
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
