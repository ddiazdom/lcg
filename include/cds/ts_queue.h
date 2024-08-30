//
// Created by Diaz, Diego on 31.3.2023.
//

#ifndef PARALLEL_PARSING_THREAD_SAFE_QUEUE_H
#define PARALLEL_PARSING_THREAD_SAFE_QUEUE_H

#include<mutex>
#include<queue>
#include<condition_variable>
#include<thread>
#include <atomic>

template<typename T>
class ts_queue {
public:
    ts_queue() = default;
    ~ts_queue() {
        done();
    };

    void push(const T& item) {
        while(q_lock.test_and_set(std::memory_order_acquire));
        m_queue.push(item);
        q_lock.clear(std::memory_order_release);
        //m_condition.notify_one();
    }

    void push(T&& item) {
        while(q_lock.test_and_set(std::memory_order_acquire));
        m_queue.push(std::move(item));
        q_lock.clear(std::memory_order_release);
        //m_condition.notify_one();
    }

    bool pop(T& item) {
        while(true){
            while(q_lock.test_and_set(std::memory_order_acquire));
            if(!m_queue.empty()){
                item = std::move(m_queue.front());
                m_queue.pop();
                q_lock.clear(std::memory_order_release);
                break;
            } else {
                if(m_done.load(std::memory_order_acquire)){
                    q_lock.clear(std::memory_order_release);
                    return false;
                }else{
                    q_lock.clear(std::memory_order_release);
                    //std::this_thread::yield();
                }
            }
        }
        return true;
        //std::unique_lock guard(q_mtx);
        //m_condition.wait(guard, [&]() { return !m_queue.empty() || m_done; });
        //if(m_done){
        //    return false;
        //}
        //item = std::move(m_queue.front());
        //m_queue.pop();
        return true;
    }

    std::size_t size() {
        size_t sz;
        while(q_lock.test_and_set(std::memory_order_acquire));
        //std::unique_lock guard(m_queue_lock);
        sz = m_queue.size();
        q_lock.clear(std::memory_order_release);
        return sz;
    }

    bool empty() {
        bool ept;
        while(q_lock.test_and_set(std::memory_order_acquire));
        //std::unique_lock guard(m_queue_lock);
        ept = m_queue.empty();
        q_lock.clear(std::memory_order_release);
        return ept;
    }

    void done() {
        //while(q_lock.test_and_set(std::memory_order_acquire)){}
        m_done.store(true, std::memory_order_release);
        //q_lock.clear(std::memory_order_release);
        //{
        //    std::unique_lock guard(m_queue_lock);
        //m_done = true;
        //}
        //m_condition.notify_all();
    }

private:
    using queue_t = std::queue<T>;
    queue_t m_queue;
    //mutable std::mutex q_mtx;
    std::atomic_flag q_lock = ATOMIC_FLAG_INIT;
    //std::condition_variable m_condition;
    std::atomic_bool m_done{false};
};
#endif //PARALLEL_PARSING_THREAD_SAFE_QUEUE_H
