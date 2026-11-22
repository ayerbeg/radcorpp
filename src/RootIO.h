// ROOT input/output of a RunInput (format "root"). Without ROOT both functions throw.
#pragma once
#include <istream>
#include <string>

#include "radcor/Input.h"

namespace radcor {

// CSV parser on a stream (defined in Input.cpp); `source` names it in messages.
RunInput parseCsv(std::istream& is, const std::string& source, const std::string& kind);

// Read TTree "points" and TNamed "run" from `file`. Branches are read by name with any numeric
// type (Double_t, Float_t, Int_t...); required: E0, nu, xs; optional as in the CSV format.
RunInput readRootInput(const std::string& file, const std::string& kind);

// Write `in` as TTree "points" (Double_t branches, Int_t config, std::string acc_table and
// title) and TNamed "run".
void writeRootInput(const RunInput& in, const std::string& file);

}  // namespace radcor
