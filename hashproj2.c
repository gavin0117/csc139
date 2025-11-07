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
   ======================================================================= */

static void *heap_start = NULL;
static node_t *free_list = NULL;
static node_t *last_search = NULL;
static int heap_initialized = 0;

/* Round size to 8-byte alignment */
size_t round_up_size(size_t size) {
    if (size % 8 == 0) {
        return size;
    }
    return size + (8 - (size % 8));
}

/* Initialize heap with one large free block */
void initialize_heap(void) {
    heap_start = init_umem();

    if (heap_start == NULL) {
        fprintf(stderr, "Failed to initialize heap\n");
        exit(1);
    }

    free_list = (node_t *)heap_start;
    free_list->size = UMEM_SIZE - sizeof(node_t);
    free_list->next = NULL;

    last_search = NULL;
    heap_initialized = 1;
}

/* Allocate memory using next-fit strategy */
void *umalloc(size_t size) {
    if (!heap_initialized) {
        initialize_heap();
    }

    if (size == 0) {
        return NULL;
    }

    size_t aligned_size = round_up_size(size);

    node_t *start = last_search ? last_search : free_list;
    node_t *current = start;
    node_t *prev = NULL;

    int wrapped = 0;

    if (current != free_list) {
        node_t *p = free_list;
        while (p != NULL && p->next != current) {
            p = p->next;
        }
        prev = p;
    }

    while (current != NULL) {
        if ((size_t)current->size >= aligned_size) {
            size_t remaining = current->size - aligned_size;

            last_search = current->next;

            if (remaining >= sizeof(node_t) + 16) {
                // Split block
                header_t *alloc_header = (header_t *)current;
                alloc_header->size = aligned_size;
                alloc_header->magic = MAGIC;

                void *new_free_addr = (char *)current + sizeof(header_t) + aligned_size;
                node_t *new_free = (node_t *)new_free_addr;

                new_free->size = remaining - sizeof(node_t);
                new_free->next = current->next;

                if (prev == NULL) {
                    free_list = new_free;
                } else {
                    prev->next = new_free;
                }

                if (last_search == current->next) {
                    last_search = new_free;
                }

                return (void *)((char *)alloc_header + sizeof(header_t));

            } else {
                // Use entire block
                header_t *alloc_header = (header_t *)current;
                alloc_header->size = current->size;
                alloc_header->magic = MAGIC;

                if (prev == NULL) {
                    free_list = current->next;
                } else {
                    prev->next = current->next;
                }

                return (void *)((char *)alloc_header + sizeof(header_t));
            }
        }

        prev = current;
        current = current->next;

        if (current == NULL && !wrapped) {
            current = free_list;
            prev = NULL;
            wrapped = 1;
        }

        if (wrapped && current == start) {
            break;
        }
    }

    return NULL;
}

/* Free memory and coalesce adjacent blocks */
void ufree(void *ptr) {
    if (ptr == NULL) {
        return;
    }

    header_t *header = (header_t *)((char *)ptr - sizeof(header_t));

    if (header->magic != MAGIC) {
        fprintf(stderr, "Error: Invalid magic number. Possible corruption or double free.\n");
        exit(1);
    }

    node_t *freed_block = (node_t *)header;
    freed_block->size = header->size;
    freed_block->next = NULL;

    node_t *current = free_list;
    node_t *prev = NULL;

    while (current != NULL && (void *)current < (void *)freed_block) {
        prev = current;
        current = current->next;
    }

    if (current == freed_block) {
        fprintf(stderr, "Error: Double free detected\n");
        exit(1);
    }

    if (prev == NULL) {
        freed_block->next = free_list;
        free_list = freed_block;
    } else {
        freed_block->next = current;
        prev->next = freed_block;
    }

    // Coalesce with next
    if (freed_block->next != NULL) {
        void *end_of_freed = (char *)freed_block + sizeof(node_t) + freed_block->size;

        if (end_of_freed == (void *)freed_block->next) {
            freed_block->size += sizeof(node_t) + freed_block->next->size;
            freed_block->next = freed_block->next->next;
        }
    }

    // Coalesce with previous
    if (prev != NULL) {
        void *end_of_prev = (char *)prev + sizeof(node_t) + prev->size;

        if (end_of_prev == (void *)freed_block) {
            prev->size += sizeof(node_t) + freed_block->size;
            prev->next = freed_block->next;
        }
    }
}


unsigned long process_block(const unsigned char *buf, size_t len) {
    unsigned long freq[SYMBOLS] = {0};

    for (size_t i = 0; i < len; i++) {
        unsigned char symbol = buf[i];
        freq[symbol]++;
    }

    Node *root = build_tree(freq);
    unsigned long h = hash_tree(root, 0);
    free_tree(root);

    return h;
}

int run_single(const char *filename) {
    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    unsigned char buffer[BLOCK_SIZE];
    unsigned long final_hash = 0;
    int block_num = 0;

    size_t bytes_read;
    while ((bytes_read = fread(buffer, 1, BLOCK_SIZE, file)) > 0) {
        unsigned long block_hash = process_block(buffer, bytes_read);
        print_intermediate(block_num, block_hash, getpid());
        final_hash = (final_hash + block_hash) % LARGE_PRIME;
        block_num++;
    }

    fclose(file);
    print_final(final_hash);

    return 0;
}

int run_multi(const char *filename) {
    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    unsigned char buffer[BLOCK_SIZE];
    unsigned long final_hash = 0;
    int block_num = 0;

    pid_t child_pids[1000];
    int num_children = 0;

    size_t bytes_read;
    while ((bytes_read = fread(buffer, 1, BLOCK_SIZE, file)) > 0) {
        int pipefd[2];
        if (pipe(pipefd) == -1) {
            perror("pipe");
            fclose(file);
            return 1;
        }

        pid_t pid = fork();

        if (pid < 0) {
            perror("fork");
            fclose(file);
            return 1;
        }

        if (pid == 0) {
            // Child process
            fclose(file);
            close(pipefd[0]);

            unsigned long block_hash = process_block(buffer, bytes_read);
            write(pipefd[1], &block_hash, sizeof(block_hash));

            close(pipefd[1]);
            _exit(0);

        } else {
            // Parent process
            close(pipefd[1]);

            unsigned long block_hash;
            ssize_t bytes = read(pipefd[0], &block_hash, sizeof(block_hash));

            if (bytes != sizeof(block_hash)) {
                fprintf(stderr, "Error reading from pipe\n");
                close(pipefd[0]);
                fclose(file);
                return 1;
            }

            close(pipefd[0]);
            print_intermediate(block_num, block_hash, pid);
            final_hash = (final_hash + block_hash) % LARGE_PRIME;

            child_pids[num_children] = pid;
            num_children++;

            block_num++;
        }
    }

    fclose(file);

    // Wait for all children
    for (int i = 0; i < num_children; i++) {
        int status;
        waitpid(child_pids[i], &status, 0);
    }

    print_final(final_hash);

    return 0;
}
