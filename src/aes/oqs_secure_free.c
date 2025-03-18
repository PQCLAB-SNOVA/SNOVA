/**
 * Add missing function needed by the oqs derived aes
 * 
 * Adapted from
 * 
 * https://github.com/open-quantum-safe/liboqs
 * commit 3488f0a598c64b730ee2e2a4acb38e1a51797c99
 * 
 * MIT license
 */

#include <stddef.h>
#include <string.h>
#include <stdlib.h>

void OQS_MEM_secure_free(void *ptr, size_t len)
{
    if (ptr != NULL)
    {
        memset(ptr, 0, len);
        free(ptr); // IGNORE free-check
    }
}

