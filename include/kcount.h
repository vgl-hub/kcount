#ifndef KCOUNT_H
#define KCOUNT_H

class KDB : public Kmap<KDB, UserInputKcount, uint8_t, uint32_t> { // CRTP
    UserInputKcount userInput;
public:
    KDB(UserInputKcount& userInput) : Kmap{userInput}, userInput(userInput) {
        DBextension = "kc";
    }
};

#endif /* KCOUNT_H */
