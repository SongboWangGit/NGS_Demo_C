//
// Created by sbwang on 19-3-26.
//

#ifndef NGS_DEMO_SEMAPHORE_H
#define NGS_DEMO_SEMAPHORE_H

#include <mutex>
#include <condition_variable>
#include <thread>
class Semaphore
{
public:
    Semaphore();
    explicit Semaphore(unsigned int count); //用无符号数表示信号量资源
    explicit Semaphore(Semaphore &semaphore1);
    ~Semaphore();

public:
    void wait();
    void signal();

    Semaphore& operator=(const Semaphore &x);
private:
    int m_count; //计数器必须是有符号数
    std::mutex m_mutex;
    std::condition_variable m_condition_variable;
};

#endif //NGS_DEMO_SEMAPHORE_H
Semaphore::Semaphore() {

}
Semaphore::Semaphore(unsigned int count) :m_count(count) {
}

Semaphore::Semaphore(Semaphore &semaphore1) {

    m_count = semaphore1.m_count;
}

Semaphore& Semaphore::operator=(const Semaphore &x) {
    this->m_count = x.m_count;
}
Semaphore::~Semaphore()
{
}

void Semaphore::wait() {
    std::unique_lock<std::mutex> unique_lock(m_mutex);
    --m_count;
    while (m_count < 0) {
        m_condition_variable.wait(unique_lock);
    }
}

void Semaphore::signal() {
    std::lock_guard<std::mutex> lg(m_mutex);
    if (++m_count < 1) {
        m_condition_variable.notify_one();
    }
}
