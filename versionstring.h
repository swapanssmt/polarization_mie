#include "version.h"
#include <stdio.h>
#include <string.h>

void version_string(char *buf) {
  
  sprintf(buf,"--------------------------------------------\n");
  sprintf(buf+strlen(buf),"  Version:  %s\n", build_version);
  sprintf(buf+strlen(buf),"  Revision: %s\n", build_revision);
  sprintf(buf+strlen(buf),"--------------------------------------------\n");
}

