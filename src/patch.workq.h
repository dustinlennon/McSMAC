22a23,24
> #ifndef BUTENHOF_H_
> #define BUTENHOF_H_
30a33
>     void                        *res;
46c49
<     void                (*engine)(void *arg);   /* user engine */
---
>     void                (*engine)(void *res, void *data);   /* user engine */
57c60,61
<     void        (*engine)(void *));     /* engine routine */
---
>     void        (*engine)(void *, void*));     /* engine routine */
> 
59c63,64
< extern int workq_add (workq_t *wq, void *data);
---
> extern int workq_add (workq_t *wq, void* res, void *data);
> #endif
