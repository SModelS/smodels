#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/LogBicubicInterpolator.h"
#include <memory>
using namespace std;

int main() {

  // Get a PDF set object and use it to get a vector of heap-allocated PDFs
  LHAPDF::PDFSet set("CT10nlo");
  const auto pdfs = set.mkPDFs<shared_ptr<const LHAPDF::PDF>>();

  // Test the built-in cache
  cout << "Original cache sizes = "
       << LHAPDF::LogBicubicInterpolator::XCaches::SIZE << " + "
       << LHAPDF::LogBicubicInterpolator::Q2Caches::SIZE << endl;
  cout << pdfs[0]->xfxQ(21, 1e-3, 126.0) << endl;

  // Try several different cache sizes at runtime
  for (size_t xcachesize : {1,2,4,8,16,32,64}) {
    for (size_t q2cachesize : {1,2,4}) {
      LHAPDF::LogBicubicInterpolator::XCaches::setup(xcachesize);
      LHAPDF::LogBicubicInterpolator::Q2Caches::setup(q2cachesize);
      cout << "New cache sizes = "
           << LHAPDF::LogBicubicInterpolator::XCaches::SIZE << " + "
           << LHAPDF::LogBicubicInterpolator::Q2Caches::SIZE << endl;
      cout << pdfs[0]->xfxQ(21, 1e-3, 126.0) << endl;

      for (double q = 10; q < 1e4; q *= 9) {
        for (double x = 1e-7; x < 1; x *= 10) {
          for (int pid = -5; pid < 6; ++pid) {
            for (auto p : pdfs) {
              const double xf_g = p->xfxQ(pid, x, q);
              cout << pid << ", " << x << ", " << q << ": " << xf_g << endl;
            }
          }
        }
      }

    }
  }

  return 0;
}
