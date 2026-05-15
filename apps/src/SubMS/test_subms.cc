#include <filesystem>
#include <SubMS/subms.h>
#include <gtest/gtest.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/casa/Arrays/Array.h>

using namespace std;
using namespace std::filesystem;
using namespace casacore;

namespace test {

const path goldDir = path(SUBMS_TEST_DATA_DIR) / "gold_standard";

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static void copyGoldMS(const path& testDir)
{
    std::filesystem::copy(goldDir / "CYGTST.corespiral.ms",
                          testDir / "test.ms",
                          copy_options::recursive);
}

// ---------------------------------------------------------------------------
// UIFactory
// ---------------------------------------------------------------------------

TEST(SubMSTest, UIFactory) {
    int argc = 4;
    char arg0[] = "./subms";
    char arg1[] = "help=noprompt";
    char arg2[] = "ms=testin.ms";
    char arg3[] = "outms=testout.ms";
    char* argv[] = {arg0, arg1, arg2, arg3};

    string MSNBuf="testin.ms", OutMSBuf="testout.ms", WhichColStr="data",
           fieldStr="*", timeStr, spwStr="*", baselineStr, uvdistStr,
           taqlStr, scanStr, arrayStr, corrStr, intentStr, obsStr, combineStr;
    double integ    = -1.0;
    int    chanStep = 1;
    bool   overwrite = false, verbose = false;
    Bool   restartUI = false, interactive = false;

    UI(restartUI, argc, argv, interactive,
       MSNBuf, OutMSBuf, WhichColStr,
       fieldStr, timeStr, spwStr, baselineStr,
       scanStr, arrayStr, uvdistStr, taqlStr,
       corrStr, intentStr, obsStr, combineStr,
       integ, chanStep, overwrite, verbose);
}

// ---------------------------------------------------------------------------
// Basic split -- output MS exists and has rows
// ---------------------------------------------------------------------------

TEST(SubMSTest, BasicSplitProducesOutputMS) {
    string testName = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    path testDir = path(SUBMS_TEST_DATA_DIR) / testName;
    std::filesystem::remove_all(testDir);
    std::filesystem::create_directory(testDir);
    copyGoldMS(testDir);
    std::filesystem::current_path(testDir);

    SubMS_func("test.ms", "out.ms");

    ASSERT_TRUE(std::filesystem::exists("out.ms"));
    MeasurementSet ms("out.ms", Table::Old);
    EXPECT_GT(ms.nrow(), 0u);

    std::filesystem::current_path(testDir.parent_path());
    std::filesystem::remove_all(testDir);
}

// ---------------------------------------------------------------------------
// Overwrite=false throws when output already exists
// ---------------------------------------------------------------------------

TEST(SubMSTest, OverwriteFalseThrowsIfOutputExists) {
    string testName = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    path testDir = path(SUBMS_TEST_DATA_DIR) / testName;
    std::filesystem::remove_all(testDir);
    std::filesystem::create_directory(testDir);
    copyGoldMS(testDir);
    std::filesystem::current_path(testDir);

    // Create a dummy output so the guard triggers
    std::filesystem::create_directory("out.ms");

    EXPECT_THROW(
        SubMS_func("test.ms", "out.ms", "data",
                   "*", "", "*", "", "", "", "", "", "", "", "", "",
                   -1.0, 1, /*overwrite=*/false),
        AipsError);

    std::filesystem::current_path(testDir.parent_path());
    std::filesystem::remove_all(testDir);
}

// ---------------------------------------------------------------------------
// Overwrite=true replaces existing output
// ---------------------------------------------------------------------------

TEST(SubMSTest, OverwriteTrueReplacesExistingOutput) {
    string testName = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    path testDir = path(SUBMS_TEST_DATA_DIR) / testName;
    std::filesystem::remove_all(testDir);
    std::filesystem::create_directory(testDir);
    copyGoldMS(testDir);
    std::filesystem::current_path(testDir);

    std::filesystem::create_directory("out.ms");

    EXPECT_NO_THROW(
        SubMS_func("test.ms", "out.ms", "data",
                   "*", "", "*", "", "", "", "", "", "", "", "", "",
                   -1.0, 1, /*overwrite=*/true));

    ASSERT_TRUE(std::filesystem::exists("out.ms"));
    MeasurementSet ms("out.ms", Table::Old);
    EXPECT_GT(ms.nrow(), 0u);

    std::filesystem::current_path(testDir.parent_path());
    std::filesystem::remove_all(testDir);
}

// ---------------------------------------------------------------------------
// chanStep=1 and chanStep=2 produce different channel counts
// ---------------------------------------------------------------------------

TEST(SubMSTest, ChanStepReducesChannelCount) {
    string testName = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    path testDir = path(SUBMS_TEST_DATA_DIR) / testName;
    std::filesystem::remove_all(testDir);
    std::filesystem::create_directory(testDir);
    copyGoldMS(testDir);
    std::filesystem::current_path(testDir);

    SubMS_func("test.ms", "out_step1.ms", "data",
               "*", "", "*", "", "", "", "", "", "", "", "", "",
               -1.0, /*chanStep=*/1);

    SubMS_func("test.ms", "out_step2.ms", "data",
               "*", "", "*", "", "", "", "", "", "", "", "", "",
               -1.0, /*chanStep=*/2);

    MeasurementSet ms1("out_step1.ms", Table::Old);
    MeasurementSet ms2("out_step2.ms", Table::Old);

    ArrayColumn<Complex> col1(ms1, MS::columnName(MS::DATA));
    ArrayColumn<Complex> col2(ms2, MS::columnName(MS::DATA));

    IPosition shape1 = col1.shape(0);
    IPosition shape2 = col2.shape(0);

    if (shape1(1) > 1) {
        // step=2 averages pairs of channels: output has half as many
        EXPECT_EQ(shape2(1), shape1(1) / 2);
    } else {
        // MS has only 1 channel -- chanStep=2 cannot reduce below 1
        EXPECT_EQ(shape2(1), 1);
    }

    std::filesystem::current_path(testDir.parent_path());
    std::filesystem::remove_all(testDir);
}

// ---------------------------------------------------------------------------
// chanStep < 1 throws before touching the MS
// ---------------------------------------------------------------------------

TEST(SubMSTest, InvalidChanStepThrows) {
    string testName = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    path testDir = path(SUBMS_TEST_DATA_DIR) / testName;
    std::filesystem::remove_all(testDir);
    std::filesystem::create_directory(testDir);
    copyGoldMS(testDir);
    std::filesystem::current_path(testDir);

    EXPECT_THROW(
        SubMS_func("test.ms", "out.ms", "data",
                   "*", "", "*", "", "", "", "", "", "", "", "", "",
                   -1.0, /*chanStep=*/0),
        AipsError);

    std::filesystem::current_path(testDir.parent_path());
    std::filesystem::remove_all(testDir);
}

// ---------------------------------------------------------------------------
// verbose=true runs without error (output goes to LogIO, not testable by value)
// ---------------------------------------------------------------------------

TEST(SubMSTest, VerboseRunsWithoutError) {
    string testName = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    path testDir = path(SUBMS_TEST_DATA_DIR) / testName;
    std::filesystem::remove_all(testDir);
    std::filesystem::create_directory(testDir);
    copyGoldMS(testDir);
    std::filesystem::current_path(testDir);

    EXPECT_NO_THROW(
        SubMS_func("test.ms", "out.ms", "data",
                   "*", "", "*", "", "", "", "", "", "", "", "", "",
                   -1.0, 1, false, /*verbose=*/true));

    EXPECT_TRUE(std::filesystem::exists("out.ms"));

    std::filesystem::current_path(testDir.parent_path());
    std::filesystem::remove_all(testDir);
}

}; // namespace test
