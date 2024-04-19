//
// Created by Diaz, Diego on 23.4.2023.
//

#ifndef SIMPLE_EXAMPLE_TS_PRIORITY_QUEUE_H
#define SIMPLE_EXAMPLE_TS_PRIORITY_QUEUE_H

#include <queue>

template<typename value_type,
         typename container_type,
         typename comparator_type>
class ts_priority_queue{

private:
    std::atomic_flag q_lock = ATOMIC_FLAG_INIT;
    std::priority_queue<value_type, container_type, comparator_type> p_queue;
    std::atomic_bool m_done{false};

public:
    ~ts_priority_queue() {
        done();
    };

    void push(const value_type& item) {
        while(q_lock.test_and_set(std::memory_order_acquire));
        p_queue.push(item);
        q_lock.clear(std::memory_order_release);
    }

    void push(value_type&& item) {
        while(q_lock.test_and_set(std::memory_order_acquire));
        p_queue.push(std::move(item));
        q_lock.clear(std::memory_order_release);
    }

    bool pop(value_type& item) {
        while(true){
            while(q_lock.test_and_set(std::memory_order_acquire));
            if(!p_queue.empty()){
                item = std::move(p_queue.top());
                p_queue.pop();
                q_lock.clear(std::memory_order_release);
                break;
            } else {
                if(m_done.load(std::memory_order_acquire)){
                    q_lock.clear(std::memory_order_release);
                    return false;
                }else{
                    q_lock.clear(std::memory_order_release);
                    std::this_thread::yield();
                }
            }
        }
        return true;
    }

    bool top(value_type& item){
        while(true){
            while(q_lock.test_and_set(std::memory_order_acquire));
            if(!p_queue.empty()){
                item = p_queue.top();
                q_lock.clear(std::memory_order_release);
                break;
            } else {
                if(m_done.load(std::memory_order_acquire)){
                    q_lock.clear(std::memory_order_release);
                    return false;
                }else{
                    q_lock.clear(std::memory_order_release);
                    std::this_thread::yield();
                }
            }
        }
        return true;
    }

    std::size_t size() {
        size_t sz;
        while(q_lock.test_and_set(std::memory_order_acquire));
        sz = p_queue.size();
        q_lock.clear(std::memory_order_release);
        return sz;
    }

    bool empty() {
        bool ept;
        while(q_lock.test_and_set(std::memory_order_acquire));
        ept = p_queue.empty();
        q_lock.clear(std::memory_order_release);
        return ept;
    }

    void done() {
        m_done.store(true, std::memory_order_release);
    }
};
#endif //SIMPLE_EXAMPLE_TS_PRIORITY_QUEUE_H
