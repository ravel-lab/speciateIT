  fprintf(stderr, "Testing existence of the limbolic link file %s\n", symblink.c_str());
  if ( exists( symblink.c_str() ) )
  {
    fprintf(stderr, "\tDetected the limbolic link file\n");
    remove(symblink.c_str());
  } else {
    fprintf(stderr, "\tThe limbolic link file not detected\n");
  }

  if ( access(symblink.c_str(), F_OK) != -1) {
    fprintf(stderr, "\tDetected the limbolic link file using access()\n");
  } else {
    fprintf(stderr, "\tThe limbolic link file not detected with access()\n");
  }
