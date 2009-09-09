###############################################
# cmake configuration file for MarlinPandora
# @author Jan Engels, DESY
###############################################

SET( MarlinPandora_FOUND FALSE )
MARK_AS_ADVANCED( MarlinPandora_FOUND )

# do not store find results in cache
SET( MarlinPandora_INCLUDE_DIR MarlinPandora_INCLUDE_DIR-NOTFOUND )

FIND_PATH( MarlinPandora_INCLUDE_DIR
	NAMES Pandora.h
	PATHS /usera/marshall/ilcsoft/MarlinPandora
	PATH_SUFFIXES include
	NO_DEFAULT_PATH
)
IF( NOT MarlinPandora_INCLUDE_DIR )
    MESSAGE( STATUS "Check for MarlinPandora: ${MarlinPandora_HOME}"
					" -- failed to find MarlinPandora include directory!!" )
ELSE( NOT MarlinPandora_INCLUDE_DIR )
    MARK_AS_ADVANCED( MarlinPandora_INCLUDE_DIR )
ENDIF( NOT MarlinPandora_INCLUDE_DIR )


# do not store find results in cache
SET( MarlinPandora_LIB MarlinPandora_LIB-NOTFOUND )

FIND_LIBRARY( MarlinPandora_LIB
	NAMES MarlinPandora
	PATHS /usera/marshall/ilcsoft/MarlinPandora
	PATH_SUFFIXES lib
	NO_DEFAULT_PATH
)
IF( NOT MarlinPandora_LIB )
    MESSAGE( STATUS "Check for MarlinPandora: ${MarlinPandora_HOME}"
					" -- failed to find MarlinPandora library!!" )
ELSE( NOT MarlinPandora_LIB )
    MARK_AS_ADVANCED( MarlinPandora_LIB )
ENDIF( NOT MarlinPandora_LIB )


# set variables and display results
IF( MarlinPandora_INCLUDE_DIR AND MarlinPandora_LIB )
    SET( MarlinPandora_FOUND TRUE )
    SET( MarlinPandora_INCLUDE_DIRS ${MarlinPandora_INCLUDE_DIR} )
    SET( MarlinPandora_LIBRARY_DIRS "/usera/marshall/ilcsoft/MarlinPandora/lib" )
	SET( MarlinPandora_LIBRARIES ${MarlinPandora_LIB} )
    MARK_AS_ADVANCED( MarlinPandora_INCLUDE_DIRS MarlinPandora_LIBRARY_DIRS MarlinPandora_LIBRARIES )
	MESSAGE( STATUS "Check for MarlinPandora: ${MarlinPandora_HOME} -- works" )
ELSE( MarlinPandora_INCLUDE_DIR AND MarlinPandora_LIB )
	IF( MarlinPandora_FIND_REQUIRED )
		MESSAGE( FATAL_ERROR "Check for MarlinPandora: ${MarlinPandora_HOME} -- failed!!" )
    ELSE( MarlinPandora_FIND_REQUIRED )
        MESSAGE( STATUS "Check for MarlinPandora: ${MarlinPandora_HOME}"
						" -- failed!! will skip this package..." )
    ENDIF( MarlinPandora_FIND_REQUIRED )
ENDIF( MarlinPandora_INCLUDE_DIR AND MarlinPandora_LIB )
