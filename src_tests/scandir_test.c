//
// This is a modification of ~/devel/speciateIT/dirent-1.23.2/tests/t-scandir.c
//

// compile with
// clang scandir_test.c -o scandir_test

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _MSC_VER
#   include <direct.h>
#else
#   include <unistd.h>
#endif
#include <sys/stat.h>
#include <dirent.h>

//#undef NDEBUG
//#include <assert.h>


/* Filter and sort functions */
static int no_directories (const struct dirent *entry);
static int reverse_alpha (const struct dirent **a, const struct dirent **b);

int main(int argc, char *argv[])
{

    if (argc != 2) {
      fprintf(stderr, "Usage: %s <dir name>\n", argv[0]);
      exit(EXIT_FAILURE);
    }

    char *dir = argv[1];

    struct dirent **files;
    // Read directory entries in alphabetic order
    int n = scandir(dir, &files, NULL, alphasort); // reverse_alpha for reverse alphasort

    fprintf(stderr, "\n\n========================================\n");
    fprintf(stderr, "Numer of files found (including subdirectories) %d\n", n);
    fprintf(stderr, "Files found\n");
    char *s;
    for ( int i = 0; i < n; i++ )
    {
      s = files[i]->d_name;
      if ( strcmp (s, ".") != 0 && strcmp (s, "..") != 0  )
      {
        if ( files[i]->d_type != DT_DIR )
          fprintf(stderr, "%s\n", s);
        else
          fprintf(stderr, "%s/\n", s);
      }
      free (files[i]);
    }

    //
    // The same, but now does not list directories and sort in reverse order
    //
    n = scandir(dir, &files, no_directories, reverse_alpha);

    fprintf(stderr, "\n\n========================================\n");
    fprintf(stderr, "Another scan excluding subdirectories\n");
    fprintf(stderr, "Numer of files found %d\n", n);
    fprintf(stderr, "Files found\n");
    for ( int i = 0; i < n; i++ )
    {
      s = files[i]->d_name;
      if ( strcmp (s, ".") != 0 && strcmp (s, "..") != 0  )
        fprintf(stderr, "%s\n", s);
      free (files[i]);
    }

    return EXIT_SUCCESS;
}

/* Filter out directories */
static int no_directories (const struct dirent *entry)
{
    int pass;

    if (entry->d_type != DT_DIR) {
        pass = 1;
    } else {
        pass = 0;
    }

    return pass;
}

/* Sort in reverse direction */
static int reverse_alpha(const struct dirent **a, const struct dirent **b)
{
    return strcoll ((*b)->d_name, (*a)->d_name);
}
