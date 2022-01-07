#!/bin/bash
# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# SUNDIALS testing driver script
# ------------------------------------------------------------------------------

# Exit when any command fails
set -e

# Help message
help ()
{
    cat <<EOF

    $0 [--help] [--buildjobs NUM] [--testjobs NUM] [--testtype TYPE]
       [--branchname NAME] [--changebranch NAME] [--tarball PACKAGE]
       [--realtype TYPE] [--indexsize SIZE] [--libtype TYPE] [--tpls]
       [--suntesttype TYPE] [--phase PHASE] [--env ENVFILE] [EXTRA INPUTS]

    where:

    --help               Displays this message and exit.
    --buildjobs NUM      Set the number of build jobs (default 1).
    --testjobs NUM       Set the number of test jobs (default 1).
    --testtype TYPE      The test type to run (default custom):
                           branch  -- quick tests with a few configurations
                           pr      -- create tarball, test more configurations
                           release -- create tarball, test more configurations
                           custom  -- single user defined test configuration
    --branchname NAME    Branch name or PR name (overrides --testtype,
                         default NONE).
    --changebranch NAME  Source branch name in a PR (overrides --testtype,
                         default NONE).
    --tarball PACKAGE    Tarball to create in a custom test (default NONE):
                           [sundials|arkode|cvode|cvodes|ida|idas|kinsol|all]
    --realtype TYPE      Real type to use in a custom test (default double):
                           [single|double|extended]
    --indexsize SIZE     Index size to use in a custom test (default 64):
                           [32|64]
    --libtype TYPE       Library type to create a custom test (default both):
                           [static|shared|both]
    --tpls               Enable TPLs in a custom test.
    --suntesttype TYPE   SUNDIALS test type for a custom test (default DEV):
                           NONE -- compile test only
                           STD  -- run standard tests
                           DEV  -- run development tests
    --phase PHASE        Testing phase to stop after:
                           ENV -- after setting up the environment
                           CONFIG -- after configuring with CMake
                           BUILD -- after make
                           TEST -- after make test
                           INSTALL -- after make install
                           TEST_INSTALL -- after make test_install
                           TEST_INSTALL_ALL -- after make test_install_all
    --env ENVFILE        Environment file to use
    EXTRA INPUTS         Extra inputs passed to the environment file

    Example usage:

    $0
    $0 --testtype release --buildjobs 4
    $0 --phase CONFIG --indexsize 32 --tpls --env env/default.sh

EOF
}

# Remove old logs
\rm -rf ./*.log

# Print input args
echo "./test_driver.sh $*" | tee -a suntest.log

# --------------
# Default values
# --------------

buildjobs=0
testjobs=0
testtype="CUSTOM"
branchname=""
changebranch=""
tarball="NONE"
realtype="double"
indexsize="64"
libtype="both"
tpls="OFF"
suntesttype="DEV"
phase=""

# ------------
# Parse inputs
# ------------

# Store extra inputs
EXTRA_ARGS=()

# Iterate over inputs
while [[ $# -gt 0 ]]; do

    input="$1"

    case $input in

        --help)
            help
            exit 1;;

        --testtype)
            testtype=$2
            case "$testtype" in
                BRANCH|branch|Branch)
                    testtype=BRANCH
                    ;;
                PR|pr)
                    testtype=PR
                    ;;
                RELEASE|release|Release)
                    testtype=RELEASE
                    ;;
                CUSTOM|custom|Custom)
                    testtype=CUSTOM
                    ;;
                *)
                    echo "ERROR: Invaid test type $testtype"
                    help
                    exit 1;;
            esac
            shift 2;;

        --buildjobs)
            buildjobs=$2
            if [[ $buildjobs -lt 1 ]]; then
                echo "ERROR: Build jobs < 1"
                help
                exit 1
            fi
            shift 2;;

        --testjobs)
            testjobs=$2
            if [[ $testjobs -lt 1 ]]; then
                echo "ERROR: Test jobs < 1"
                help
                exit 1
            fi
            shift 2;;

        --branchname)
            branchname=$2
            shift 2;;

        --changebranch)
            changebranch=$2
            shift 2;;

        --tarball)
            tarball=$2
            case "$tarball" in
                sundials|arkode|cvode|cvodes|ida|idas|kinsol|all)
                ;;
                *)
                    echo "ERROR: Invaid tarball option $tarball"
                    help
                    exit 1;;
            esac
            shift 2;;

        --realtype)
            realtype=$2
            case "$realtype" in
                SINGLE|Single|single)
                    realtype=single
                    ;;
                DOUBLE|Double|double)
                    realtype=double
                    ;;
                EXTENDED|Extended|extended)
                    realtype=extended
                    ;;
                *)
                    echo "ERROR: Invaid real type option $realtype"
                    help
                    exit 1;;
            esac
            shift 2;;

        --indexsize)
            indexsize=$2
            case "$indexsize" in
                32|64)
                ;;
                *)
                    echo "ERROR: Invaid index size option $indexsize"
                    help
                    exit 1;;
            esac
            shift 2;;

        --libtype)
            libtype=$2
            case "$libtype" in
                STATIC|Static|static)
                    libtype=static
                    ;;
                SHARED|Shared|shared)
                    libtype=shared
                    ;;
                BOTH|Both|both)
                    libtype=both
                    ;;
                *)
                    echo "ERROR: Invaid library type option $libtype"
                    help
                    exit 1;;
            esac
            shift 2;;

        --tpls)
            tpls="ON"
            shift;;

        --suntesttype)
            suntesttype=$2
            case "$suntesttype" in
                DEV|Dev|dev)
                    suntesttype=DEV
                    ;;
                STD|Std|std)
                    suntesttype=STD
                    ;;
                NONE|None|none)
                    suntesttype=NONE
                    ;;
                *)
                    echo "ERROR: Invaid SUNDIALS test type option $suntesttype"
                    help
                    exit 1;;
            esac
            shift 2;;

        --phase)
            phase=$2
            case "$phase" in
                ENV|Env|env)
                    phase=ENV
                    ;;
                CONFIG|Config|config)
                    phase=CONFIG
                    ;;
                BUILD|Build|build)
                    phase=BUILD
                    ;;
                TEST|Test)
                    phase=TEST
                    ;;
                INSTALL|Install|install)
                    phase=INSTALL
                    ;;
                TEST_INSTALL|Test_Install|test_install)
                    phase=TEST_INSTALL
                    ;;
                TEST_INSTALL_ALL|Test_Install_All|test_install_all)
                    phase=TEST_INSTALL_ALL
                    ;;
                *)
                    echo "ERROR: Invaid phase option $phase"
                    help
                    exit 1;;
            esac
            shift 2;;

        --env)
            envfile=$2
            if [[ ! -f "$envfile" ]]; then
                echo "ERROR: Environment file does not exist $envfile"
                help
                exit 1
            else
                export SUNDIALS_ENV_FILE="${envfile}"
            fi
            shift 2;;

        *)
            EXTRA_ARGS+=("$input")
            shift;;
    esac
done

# Reset positional inputs to extra inputs
set -- "${EXTRA_ARGS[@]}"

# Update the test type based on the branch name or pull-request ID (PR-#). The
# change branch check below will update the test type for a release branch PR.
if [ -n "$branchname" ]; then
    if [[ "${branchname:0:2}" == "PR" ]]; then
        testtype="PR"
    elif [[ "${branchname:0:7}" == "RELEASE" ||
            "${branchname:0:7}" == "Release" ||
            "${branchname:0:7}" == "release" ]]; then
        testtype="RELEASE"
    else
        testtype="BRANCH"
    fi
fi

# Update the test type based on the change branch name. If this is a release
# branch, the test type is set to RELEASE to ensure a PR for a release branch
# runs the release tests rather than the PR tests.
if [ -n "$changebranch" ]; then
    if [[ "${changebranch:0:7}" == "RELEASE" ||
          "${changebranch:0:7}" == "Release" ||
          "${changebranch:0:7}" == "release" ]]; then
        testtype="RELEASE"
    fi
fi

# Initialize failure counter (0 = success)
passfail=0

# Location of testing directory
testroot=$(pwd)

# ------------------------------------------------------------------------------
# Print test header
# ------------------------------------------------------------------------------

echo "--------------------------------------------------" | tee -a suntest.log
echo "SUNDIALS $testtype test: $branchname " | tee -a suntest.log
date | tee -a suntest.log
echo "--------------------------------------------------" | tee -a suntest.log
git log -1 --pretty='Commit: %H%nAuthor: %an <%ae>%nDate:   %ad%n%n%s' | tee -a suntest.log
echo "--------------------------------------------------" | tee -a suntest.log

# ------------------------------------------------------------------------------
# Create an array of test configurations
# ------------------------------------------------------------------------------

args_realtypes=()
args_indexsizes=()
args_libtypes=()
args_tpls=()
args_suntests=()

case "$testtype" in

    BRANCH)
        # Don't creat tarballs
        tarball=NONE

        # Compile tests
        for rt in single double extended; do
            for is in 32 64; do
                args_realtypes+=("${rt}")
                args_indexsizes+=("${is}")
                args_libtypes+=("both")
                args_tpls+=("OFF")
                args_suntests+=("NONE")
            done
        done

        # Development tests
        for is in 32 64; do
            args_realtypes+=("double")
            args_indexsizes+=("${is}")
            args_libtypes+=("both")
            args_tpls+=("ON")
            args_suntests+=("DEV")
        done
        ;;

    PR)
        # Create sundials tarball
        tarball=sundials

        for rt in single double extended; do
            for is in 32 64; do
                args_realtypes+=("${rt}")
                args_indexsizes+=("${is}")
                args_libtypes+=("static")
                args_tpls+=("ON")
                # Development test output files created with double
                if [[ "${rt}" == "double" ]]; then
                    args_suntests+=("DEV")
                else
                    args_suntests+=("STD")
                fi
            done
        done
        ;;

    RELEASE)
        # Create sundials tarball
        tarball=sundials

        for rt in single double extended; do
            for is in 32 64; do
                for lt in static shared; do
                    args_realtypes+=("${rt}")
                    args_indexsizes+=("${is}")
                    args_libtypes+=("${lt}")
                    args_tpls+=("ON")
                    # Development test output files created with double
                    if [[ "${rt}" == "double" ]]; then
                        args_suntests+=("DEV")
                    else
                        args_suntests+=("STD")
                    fi
                done
            done
        done
        ;;

    CUSTOM)
        # Use default or user defined values
        args_realtypes+=("${realtype}")
        args_indexsizes+=("${indexsize}")
        args_libtypes+=("${libtype}")
        args_tpls+=("${tpls}")
        args_suntests+=("${suntesttype}")
        if [ "${realtype}" != "double" ] && [ "${suntesttype}" == "DEV" ]; then
            echo "WARNING: DEV tests may fail with ${realtype} precision"
        fi
        ;;

    *)
        echo "ERROR: Unknown test type"
        exit 1
        ;;
esac

# Number of configurations to test
nconfigs=${#args_suntests[@]}

# ------------------------------------------------------------------------------
# Create and extract tarballs (if necessary)
# ------------------------------------------------------------------------------

# Directories to run tests in (one for each tarball created)
testdir=()

if [ "$tarball" != NONE ]; then

    # Setup the environment with the first test configuration. This is really in
    # case the environment needs to be setup to build the documentation in the
    # tarscript e.g., activates a Python virtual environment.
    env_config=("${args_realtypes[0]}"
                "${args_indexsizes[0]}"
                "${args_libtypes[0]}"
                "${args_tpls[0]}"
                "${args_suntests[0]}")

    # Setup environment (creates configure.log)
    time source env/setup_env.sh "${env_config[@]}" "${EXTRA_ARGS[@]}"

    rc=${PIPESTATUS[0]}
    echo -e "\nsetup_env.sh returned $rc\n" | tee -a configure.log
    if [ "$rc" -ne 0 ]; then exit 1; fi

    # Remove old tarballs directory
    \rm -rf "$testroot/tarballs"

    # Create directory for tarballs
    mkdir "$testroot/tarballs"

    # Run tarscript to create tarballs
    cd ..
    cd scripts

    echo "START TARSCRIPT"
    time ./tarscript $tarball | tee -a tar.log

    # Check tarscript return code
    rc=${PIPESTATUS[0]}
    echo -e "\ntarscript returned $rc\n" | tee -a tar.log
    if [ "$rc" -ne 0 ]; then exit 1; fi

    # Relocate log and tarballs
    mv tar.log "$testroot/tarballs/."
    mv ../tarballs/* "$testroot/tarballs/."

    # Move to tarball directory
    cd "$testroot/tarballs"

    # Loop over tarballs and test each one
    for tb in *.tar.gz; do

        # Get package name
        package=${tb%.tar.gz}

        echo "START UNTAR"
        tar -xvzf "$tb" 2>&1 | tee -a tar.log

        # Check tar return code
        rc=${PIPESTATUS[0]}
        echo -e "\ntar -xzvf returned $rc\n" | tee -a tar.log
        if [ "$rc" -ne 0 ]; then exit 1; fi

        # Move log to package directory
        mv tar.log "$package/."

        # Copy environment and testing scripts from original test directory
        cp -r "$testroot/env" "$package/test/."
        cp "$testroot/config_cmake.py" "$package/test/."

        testdir+=("$testroot/tarballs/$package/test")
    done

else

    # No tarballs created, use the original test directory
    testdir+=("$testroot")

fi

# Number of test directories
ntestdirs=${#testdir[@]}

# ------------------------------------------------------------------------------
# Run tests
# ------------------------------------------------------------------------------

for ((j=0;j<ntestdirs;j++)); do

    cd "${testdir[j]}"

    for ((i=0;i<nconfigs;i++)); do

        # ---------------------
        # Setup the environment
        # ---------------------

        env_config=("${args_realtypes[i]}"
                    "${args_indexsizes[i]}"
                    "${args_libtypes[i]}"
                    "${args_tpls[i]}"
                    "${args_suntests[i]}")

        # Print test header for Jenkins section collapsing
        echo "TEST: ${env_config[*]}"

        # Setup environment (creates configure.log)
        time source env/setup_env.sh "${env_config[@]}" "${EXTRA_ARGS[@]}"

        rc=${PIPESTATUS[0]}
        echo -e "\nsetup_env.sh returned $rc\n" | tee -a configure.log
        if [ "$rc" -ne 0 ]; then passfail=1; break; fi

        # Check if the environment sets the number of build and test jobs but do
        # not override the command line options. If neither are set then set the
        # default value.
        if [[ $buildjobs -lt 1 ]]; then
            if [[ -n "${SUNDIALS_BUILD_JOBS}" ]]; then
                # non-empty contiguous string of digits
                re='^[0-9]+$'
                if [[ "${SUNDIALS_BUILD_JOBS}" =~ $re ]] ; then
                    buildjobs="${SUNDIALS_BUILD_JOBS}"
                else
                    buildjobs=1
                fi
            else
                buildjobs=1
            fi
        fi

        if [[ $testjobs -lt 1 ]]; then
            if [[ -n "${SUNDIALS_TEST_JOBS}" ]]; then
                # non-empty contiguous string of digits
                re='^[0-9]+$'
                if [[ "${SUNDIALS_TEST_JOBS}" =~ $re ]] ; then
                    testjobs="${SUNDIALS_TEST_JOBS}"
                else
                    testjobs=1
                fi
            else
                testjobs=1
            fi
        fi

        # Check if this is the last phase
        if [ "${phase}" == "ENV" ]; then
            echo "PASSED: ${env_config[*]}"
            continue
        fi

        # -----------------------
        # Create test directories
        # -----------------------

        tmp="${env_config[*]}"
        tmp=${tmp//" "/"_"}

        builddir="build_${tmp}"
        installdir="$(pwd)/install_${tmp}"

        # Add host name to directory names
        if [ -n "$HOST" ]; then
            builddir="${builddir}_${HOST}"
            installdir="${installdir}_${HOST}"
        elif [ -n "$HOSTNAME" ]; then
            builddir="${builddir}_${HOSTNAME}"
            installdir="${installdir}_${HOSTNAME}"
        fi

        # Remove old build and install directories
        \rm -rf "$builddir" "$installdir"

        # Create and move to new build directory, move configure log
        mkdir "$builddir"
        cd "$builddir"
        mv ../configure.log .

        # -----------------------
        # Create CMake cache file
        # -----------------------

        time python3 ../config_cmake.py \
            --readenv \
            --install-prefix "${installdir}" | tee -a configure.log

        # ---------
        # Configure
        # ---------

        echo "START CMAKE"

        time cmake -C sundials.cmake ../../. | tee -a configure.log

        rc=${PIPESTATUS[0]}
        echo -e "\ncmake returned $rc\n" | tee -a configure.log
        if [ "$rc" -ne 0 ]; then passfail=1; break; fi

        # Check if this is the last phase
        if [ "${phase}" == "CONFIG" ]; then
            echo "PASSED: ${env_config[*]}"
            cd ..
            continue
        fi

        # -----
        # Build
        # -----

        echo "START MAKE"

        time make -j "$buildjobs" 2>&1 | tee make.log

        rc=${PIPESTATUS[0]}
        echo -e "\nmake returned $rc\n" | tee -a make.log
        if [ "$rc" -ne 0 ]; then passfail=1; break; fi

        # Check if tests should be skipped (compile check only)
        if [ "${args_suntests[i]}" == "NONE" ]; then
            echo "PASSED: ${env_config[*]}"
            cd ..
            continue
        fi

        # Check if this is the last phase
        if [ "${phase}" == "BUILD" ]; then
            echo "PASSED: ${env_config[*]}"
            cd ..
            continue
        fi

        # ----
        # Test
        # ----

        echo "START TEST"

        time ctest -j "$testjobs" test 2>&1 | tee test.log

        rc=${PIPESTATUS[0]}
        echo -e "\nmake test returned $rc\n" | tee -a test.log
        if [ "$rc" -ne 0 ]; then passfail=1; break; fi

        # Check if this is the last phase
        if [ "${phase}" == "TEST" ]; then
            echo "PASSED: ${env_config[*]}"
            cd ..
            continue
        fi

        # -------
        # Install
        # -------

        echo "START INSTALL"

        time make -j "$buildjobs" install 2>&1 | tee install.log

        rc=${PIPESTATUS[0]}
        echo -e "\nmake install returned $rc\n" | tee -a install.log
        if [ "$rc" -ne 0 ]; then passfail=1; break; fi

        # Check if this is the last phase
        if [ "${phase}" == "INSTALL" ]; then
            echo "PASSED: ${env_config[*]}"
            cd ..
            continue
        fi

        # ------------
        # Test install
        # ------------

        echo "START TEST_INSTALL"

        time make test_install 2>&1 | tee test_install.log

        rc=${PIPESTATUS[0]}
        echo -e "\nmake test_install returned $rc\n" | tee -a test_install.log
        if [ "$rc" -ne 0 ]; then passfail=1; break; fi

        # Check if this is the last phase
        if [ "${phase}" == "TEST_INSTALL" ]; then
            echo "PASSED: ${env_config[*]}"
            cd ..
            continue
        fi

        echo "START TEST_INSTALL_ALL"

        time make test_install_all 2>&1 | tee test_install_all.log

        rc=${PIPESTATUS[0]}
        echo -e "\nmake test_install_all returned $rc\n" | tee -a test_install_all.log
        if [ "$rc" -ne 0 ]; then passfail=1; break; fi

        # Check if this is the last phase
        if [ "${phase}" == "TEST_INSTALL_ALL" ]; then
            echo "PASSED: ${env_config[*]}"
            cd ..
            continue
        fi

        # ------------
        # End of tests
        # ------------

        echo "PASSED: ${env_config[*]}" | tee -a suntest.log
        cd ..

    done

    if [ $passfail -ne 0 ]; then
        echo "FAILED: ${env_config[*]}" | tee -a suntest.log
        cd "${testroot}"
        break
    fi

done

# ------------------------------------------------------------------------------
# Print test footer
# ------------------------------------------------------------------------------

if [ $passfail -eq 0 ]; then
    echo "All tests PASSED" | tee -a suntest.log
fi

echo "--------------------------------------------------" | tee -a suntest.log
echo "SUNDIALS $testtype test: $branchname " | tee -a suntest.log
date | tee -a suntest.log
echo "--------------------------------------------------" | tee -a suntest.log

# ------------------------------------------------------------------------------
# Exit
# ------------------------------------------------------------------------------

if [ $passfail -eq 0 ]; then
    exit 0  # Pass
else
    exit 1  # Fail
fi
