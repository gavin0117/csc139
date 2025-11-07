#ifndef GRADING_MODE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>
#include <sys/wait.h>

#define BLOCK_SIZE 1024
#define SYMBOLS 256
#define LARGE_PRIME 2147483647   // for modular hash
#define UMEM_SIZE (128 * 1024)   // 128 KB managed heap for Step 2

/* Forward declarations for student functions*/
void *umalloc(size_t size);
void ufree(void *ptr);
unsigned long process_block(const unsigned char *buf, size_t len);
int run_single(const char *filename);
int run_multi(const char *filename);


/* =======================================================================
   PROVIDED CODE — DO NOT MODIFY
   -----------------------------------------------------------------------
   This section implements the complete Huffman tree construction logic
   and internal data structures.  It uses umalloc() and ufree() for all
   dynamic memory so that your allocator can be substituted later without
   touching this code.
   ======================================================================= */

/* ============================ STEP 2 ONLY ==============================
   The function and structures below will be used in Step 2 when you build
   your own memory allocator.  DO NOT USE OR MODIFY THIS CODE IN STEP 1.

   In Step 2, you will:
     - Use init_umem() to allocate a contiguous memory region.
     - Implement a FIRST-FIT allocator using the structures below.
     - Manage all memory inside this region manually.
   ======================================================================= */

void *init_umem(void) {
    void *ptr = mmap(NULL, UMEM_SIZE,
                     PROT_READ | PROT_WRITE,
                     MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (ptr == MAP_FAILED) {
        perror("mmap");
        return NULL;
    }
    return ptr;
}

#define MAGIC 0xDEADBEEFLL  // integrity check pattern

typedef struct {
    long size;   // Size of the block (payload only)
    long magic;  // Magic number for integrity check
} header_t;

typedef struct __node_t {
    long size;               // Size of the free block
    struct __node_t *next;   // Pointer to the next free block
} node_t;

/* =======================================================================
   Huffman and Heap Infrastructure (Given)
   ======================================================================= */

typedef struct Node {
    unsigned char symbol;
    unsigned long freq;
    struct Node *left, *right;
} Node;

typedef struct {
    Node **data;
    int size;
    int capacity;
} MinHeap;

MinHeap *heap_create(int capacity) {
    MinHeap *h = umalloc(sizeof(MinHeap));
    h->data = umalloc(sizeof(Node *) * capacity);
    h->size = 0;
    h->capacity = capacity;
    return h;
}

void heap_swap(Node **a, Node **b) {
    Node *tmp = *a; *a = *b; *b = tmp;
}

void heap_push(MinHeap *h, Node *node) {
    int i = h->size++;
    h->data[i] = node;
    while (i > 0) {
        int p = (i - 1) / 2;
        if (h->data[p]->freq < h->data[i]->freq) break;
        heap_swap(&h->data[p], &h->data[i]);
        i = p;
    }
}

Node *heap_pop(MinHeap *h) {
    if (h->size == 0) return NULL;
    Node *min = h->data[0];
    h->data[0] = h->data[--h->size];
    int i = 0;
    while (1) {
        int l = 2 * i + 1, r = l + 1, smallest = i;
        if (l < h->size && h->data[l]->freq < h->data[smallest]->freq) smallest = l;
        if (r < h->size && h->data[r]->freq < h->data[smallest]->freq) smallest = r;
        if (smallest == i) break;
        heap_swap(&h->data[i], &h->data[smallest]);
        i = smallest;
    }
    return min;
}

void heap_free(MinHeap *h) {
    ufree(h->data);
    ufree(h);
}

Node *new_node(unsigned char sym, unsigned long freq, Node *l, Node *r) {
    Node *n = umalloc(sizeof(Node));
    n->symbol = sym;
    n->freq = freq;
    n->left = l;
    n->right = r;
    return n;
}

void free_tree(Node *n) {
    if (!n) return;
    free_tree(n->left);
    free_tree(n->right);
    ufree(n);
}

Node *build_tree(unsigned long freq[SYMBOLS]) {
    MinHeap *h = heap_create(SYMBOLS);
    for (int i = 0; i < SYMBOLS; i++)
        if (freq[i] > 0)
            heap_push(h, new_node((unsigned char)i, freq[i], NULL, NULL));
    if (h->size == 0) {
        heap_free(h);
        return NULL;
    }
    while (h->size > 1) {
        Node *a = heap_pop(h);
        Node *b = heap_pop(h);
        Node *p = new_node(0, a->freq + b->freq, a, b);
        heap_push(h, p);
    }
    Node *root = heap_pop(h);
    heap_free(h);
    return root;
}

unsigned long hash_tree(Node *n, unsigned long hash) {
    if (!n) return hash;
    hash = (hash * 31 + n->freq + n->symbol) % LARGE_PRIME;
    hash = hash_tree(n->left, hash);
    hash = hash_tree(n->right, hash);
    return hash;
}

/* =======================================================================
   PROVIDED PRINTING UTILITIES (DO NOT MODIFY)
   -----------------------------------------------------------------------
   These helper functions standardize all program output for testing and
   grading. Students must call these rather than using printf() directly
   for block or final output.
   ======================================================================= */

void print_intermediate(int block_num, unsigned long hash, pid_t pid) {
#ifdef DEBUG
#  if DEBUG == 2
    printf("[PID %d] Block %d hash: %lu\n", pid, block_num, hash);
#  elif DEBUG == 1
    printf("Block %d hash: %lu\n", block_num, hash);
#  endif
#else
    (void)block_num;
    (void)hash;
    (void)pid;
#endif
}

void print_final(unsigned long final_hash) {
    printf("Final signature: %lu\n", final_hash);
}

/* =======================================================================
   MAIN DISPATCH FUNCTION
   -----------------------------------------------------------------------
   Usage:
       ./hashproj <file>          -> run single-process version
       ./hashproj <file> -m       -> run multi-process version
   ======================================================================= */

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <file> [-m]\n", argv[0]);
        return 1;
    }

    const char *filename = argv[1];
    int use_multi = (argc >= 3 && strcmp(argv[2], "-m") == 0);

    if (use_multi)
        return run_multi(filename);
    else
        return run_single(filename);
}

#endif //GRADING_MODE

/* =======================================================================
   STUDENT SECTION — STEP 2: Custom Memory Allocator
   =======================================================================
   
   Memory Layout Explanation:
   
   FREE BLOCK:
   [node_t header: 16 bytes] [payload: node_t->size bytes]
   
   ALLOCATED BLOCK:
   [header_t header: 16 bytes] [payload: header_t->size bytes]
   
   Both structures are 16 bytes, so when we allocate, we replace node_t
   with header_t in place.
   ======================================================================= */

// Pointer to the start of our managed heap
static void *heap_start = NULL;

// Pointer to the head of the free list
static node_t *free_list = NULL;

// Flag to track initialization
static int heap_initialized = 0;

/* -------------------------------------------------------------------
   Helper: Round size up to nearest 8-byte multiple for alignment
------------------------------------------------------------------- */
size_t round_up_size(size_t size) {
    if (size % 8 == 0) {
        return size;
    }
    return size + (8 - (size % 8));
}

/* -------------------------------------------------------------------
   Initialize the heap with one large free block
------------------------------------------------------------------- */
void initialize_heap(void) {
    heap_start = init_umem();
    
    if (heap_start == NULL) {
        fprintf(stderr, "Failed to initialize heap\n");
        exit(1);
    }
    
    // Create one large free block covering the entire region
    free_list = (node_t *)heap_start;
    // The size is total space minus the node_t header itself
    free_list->size = UMEM_SIZE - sizeof(node_t);
    free_list->next = NULL;
    
    heap_initialized = 1;
}

/* -------------------------------------------------------------------
   umalloc: Allocate memory using First-Fit
   
   Key insight: Both node_t and header_t are 16 bytes.
   When we allocate, we replace the node_t with a header_t.
   
   We need node_t->size >= aligned_size to have enough room.
------------------------------------------------------------------- */
void *umalloc(size_t size) {
    // Initialize on first call
    if (!heap_initialized) {
        initialize_heap();
    }
    
    // Handle zero-size requests
    if (size == 0) {
        return NULL;
    }
    
    // Round up for alignment
    size_t aligned_size = round_up_size(size);
    
    // Search the free list for a block large enough
    node_t *current = free_list;
    node_t *prev = NULL;
    
    while (current != NULL) {
        // Check if this free block has enough space
        // We need at least aligned_size bytes in the payload
        // Cast to size_t to avoid signed/unsigned comparison warning
        if ((size_t)current->size >= aligned_size) {
            // Found a suitable block!
            
            // Calculate how much space would be left after allocation
            size_t remaining = current->size - aligned_size;
            
            // Decide whether to split the block
            // Only split if remaining space can hold another node_t plus some data
            if (remaining >= sizeof(node_t) + 16) {
                // SPLIT THE BLOCK
                
                // Convert current to an allocated block
                header_t *alloc_header = (header_t *)current;
                alloc_header->size = aligned_size;
                alloc_header->magic = MAGIC;
                
                // Create a new free block from the leftover space
                // New free block starts right after our allocated payload
                void *new_free_addr = (char *)current + sizeof(header_t) + aligned_size;
                node_t *new_free = (node_t *)new_free_addr;
                
                // The new free block's payload is what's left
                new_free->size = remaining - sizeof(node_t);
                new_free->next = current->next;
                
                // Update free list
                if (prev == NULL) {
                    // We split the head of the list
                    free_list = new_free;
                } else {
                    // We split a middle/end block
                    prev->next = new_free;
                }
                
                // Return pointer to the payload
                return (void *)((char *)alloc_header + sizeof(header_t));
                
            } else {
                // DON'T SPLIT - use entire block
                
                // Convert to allocated block
                header_t *alloc_header = (header_t *)current;
                alloc_header->size = current->size;
                alloc_header->magic = MAGIC;
                
                // Remove from free list
                if (prev == NULL) {
                    // Removing head
                    free_list = current->next;
                } else {
                    // Removing from middle/end
                    prev->next = current->next;
                }
                
                // Return pointer to payload
                return (void *)((char *)alloc_header + sizeof(header_t));
            }
        }
        
        // Move to next block
        prev = current;
        current = current->next;
    }
    
    // No suitable block found
    return NULL;
}

/* -------------------------------------------------------------------
   ufree: Free memory and coalesce with adjacent blocks
------------------------------------------------------------------- */
void ufree(void *ptr) {
    // NULL pointer is valid - just do nothing
    if (ptr == NULL) {
        return;
    }
    
    // Find the header (it's right before the payload)
    header_t *header = (header_t *)((char *)ptr - sizeof(header_t));
    
    // Verify magic number
    if (header->magic != MAGIC) {
        fprintf(stderr, "Error: Invalid magic number. Possible corruption or double free.\n");
        exit(1);
    }
    
    // Convert to a free block
    node_t *freed_block = (node_t *)header;
    freed_block->size = header->size;
    freed_block->next = NULL;
    
    // Insert into free list in address order (makes coalescing easier)
    node_t *current = free_list;
    node_t *prev = NULL;
    
    // Find where to insert (keeping list sorted by address)
    while (current != NULL && (void *)current < (void *)freed_block) {
        prev = current;
        current = current->next;
    }
    
    // Check for double free
    if (current == freed_block) {
        fprintf(stderr, "Error: Double free detected\n");
        exit(1);
    }
    
    // Insert the block
    if (prev == NULL) {
        // Insert at head
        freed_block->next = free_list;
        free_list = freed_block;
    } else {
        // Insert after prev
        freed_block->next = current;
        prev->next = freed_block;
    }
    
    // Try to coalesce with next block
    if (freed_block->next != NULL) {
        // Calculate where the next block should be if they're adjacent
        void *end_of_freed = (char *)freed_block + sizeof(node_t) + freed_block->size;
        
        if (end_of_freed == (void *)freed_block->next) {
            // Adjacent! Merge them
            freed_block->size += sizeof(node_t) + freed_block->next->size;
            freed_block->next = freed_block->next->next;
        }
    }
    
    // Try to coalesce with previous block
    if (prev != NULL) {
        // Calculate where freed_block should be if prev is adjacent
        void *end_of_prev = (char *)prev + sizeof(node_t) + prev->size;
        
        if (end_of_prev == (void *)freed_block) {
            // Adjacent! Merge them
            prev->size += sizeof(node_t) + freed_block->size;
            prev->next = freed_block->next;
        }
    }
}


unsigned long process_block(const unsigned char *buf, size_t len) {
    /* Step 1: Initialize frequency table for all 256 possible symbols */
    unsigned long freq[SYMBOLS] = {0};

    /* Step 2: Count symbol frequencies in this block */
    for (size_t i = 0; i < len; i++) {
        unsigned char symbol = buf[i];
        freq[symbol]++;
    }

    /* Step 3: Build Huffman tree */
    Node *root = build_tree(freq);

    /* Step 4: Compute hash of the tree */
    unsigned long h = hash_tree(root, 0);

    /* Step 5: Free the Huffman tree */
    free_tree(root);

    /* Step 6: Return the computed hash */
    return h;
}

int run_single(const char *filename) {
    // Open file in binary mode
    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    // Buffer for reading blocks
    unsigned char buffer[BLOCK_SIZE];
    
    // Final hash accumulator
    unsigned long final_hash = 0;
    
    // Block counter
    int block_num = 0;
    
    // Read and process file in chunks
    size_t bytes_read;
    while ((bytes_read = fread(buffer, 1, BLOCK_SIZE, file)) > 0) {
        // Process this block
        unsigned long block_hash = process_block(buffer, bytes_read);
        
        // Print intermediate result
        print_intermediate(block_num, block_hash, getpid());
        
        // Accumulate into final hash
        final_hash = (final_hash + block_hash) % LARGE_PRIME;
        
        block_num++;
    }
    
    // Close file
    fclose(file);
    
    // Print final signature
    print_final(final_hash);
    
    return 0;
}

int run_multi(const char *filename) {
    // Open file in binary mode
    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    unsigned char buffer[BLOCK_SIZE];
    unsigned long final_hash = 0;
    int block_num = 0;
    
    // Array to store child PIDs so we can wait for them later
    pid_t child_pids[1000];
    int num_children = 0;
    
    size_t bytes_read;
    while ((bytes_read = fread(buffer, 1, BLOCK_SIZE, file)) > 0) {
        // Create a pipe for this child to send back its result
        int pipefd[2];
        if (pipe(pipefd) == -1) {
            perror("pipe");
            fclose(file);
            return 1;
        }
        
        // Fork a new process for this block
        pid_t pid = fork();
        
        if (pid < 0) {
            perror("fork");
            fclose(file);
            return 1;
        }
        
        if (pid == 0) {
            // CHILD PROCESS
            // Close file descriptor - child doesn't need it
            fclose(file);
            
            // Close the read end of pipe - we only write
            close(pipefd[0]);
            
            // Process this block and compute its hash
            unsigned long block_hash = process_block(buffer, bytes_read);
            
            // Write the hash to the pipe
            write(pipefd[1], &block_hash, sizeof(block_hash));
            
            // Close write end and exit immediately
            close(pipefd[1]);
            _exit(0);
            
        } else {
            // PARENT PROCESS
            // Close the write end - we only read
            close(pipefd[1]);
            
            // Read the hash from the child
            unsigned long block_hash;
            ssize_t bytes = read(pipefd[0], &block_hash, sizeof(block_hash));
            
            if (bytes != sizeof(block_hash)) {
                fprintf(stderr, "Error reading from pipe\n");
                close(pipefd[0]);
                fclose(file);
                return 1;
            }
            
            // Close read end now that we got the data
            close(pipefd[0]);
            
            // Print intermediate result with the child's PID
            print_intermediate(block_num, block_hash, pid);
            
            // Accumulate into final hash
            final_hash = (final_hash + block_hash) % LARGE_PRIME;
            
            // Save the child PID so we can wait for it
            child_pids[num_children] = pid;
            num_children++;
            
            block_num++;
        }
    }
    
    // Close the file
    fclose(file);
    
    // Wait for all child processes to finish
    for (int i = 0; i < num_children; i++) {
        int status;
        waitpid(child_pids[i], &status, 0);
    }
    
    // Print the final signature
    print_final(final_hash);
    
    return 0;
}
