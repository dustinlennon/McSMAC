27a28,31
> #ifdef __MACH__
> #include <mach/clock.h>
> #include <mach/mach.h>
> #endif
54,55c58,69
<         clock_gettime (CLOCK_REALTIME, &timeout);
<         timeout.tv_sec += 2;
---
> #ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
> clock_serv_t cclock;
> mach_timespec_t mts;
> host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
> clock_get_time(cclock, &mts);
> mach_port_deallocate(mach_task_self(), cclock);
> timeout.tv_sec = mts.tv_sec;
> timeout.tv_nsec = mts.tv_nsec;
> #else
>         clock_gettime(CLOCK_REALTIME, &timeout);
> #endif
>         timeout.tv_sec += 0.25;
103c117
<             wq->engine (we->data);
---
>             wq->engine (we->res, we->data);
153c167
< int workq_init (workq_t *wq, int threads, void (*engine)(void *arg))
---
> int workq_init (workq_t *wq, int threads, void (*engine)(void *res, void *data))
250c264
< int workq_add (workq_t *wq, void *element)
---
> int workq_add (workq_t *wq, void* res, void *element)
265a280
>     item->res = res;
