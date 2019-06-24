file(STRINGS "ppidd_SHA1" ppidd_GIT_TAG)
    FetchContent_Declare(
            ppidd_dep
            GIT_REPOSITORY https://bitbucket.org/pjknowles/ppidd
            GIT_TAG ${ppidd_GIT_TAG}
    )
