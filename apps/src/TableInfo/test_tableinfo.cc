#include <filesystem>
#include <TableInfo/tableinfo.h>
#include <gtest/gtest.h>

using namespace std;
using namespace std::filesystem;


namespace test{
TEST(TableInfoTest, UIFactory) {
    // The Factory Settings.
  int argc = 3;
  const char* argv[] = {"./tableinfo", "help=noprompt", "table=test.ms"};

  // The Factory Settings.
      bool restartUI = false;
      bool verbose = 0;
      string MSNBuf,OutBuf;
      MSNBuf="test.ms";
      OutBuf="";
      bool interactive = false;
      
      UI(restartUI,argc, (char **)argv, interactive, MSNBuf,OutBuf,verbose);
}

};
