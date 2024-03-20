# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 20:30:50 2023

@author: cosbo
"""

import telebot

import computing

from telebot.util import quick_markup


try:
    #with open("callid.txt") as f:
    #    call_id = int(f.read())

    with open("TOKEN.txt") as ff:
        API_TOKEN = ff.read()
    bot = telebot.TeleBot(API_TOKEN)
except FileNotFoundError:
    pass


chat_id = False


# Handle '/start' and '/help'
@bot.message_handler(commands=["help", "start"])
def send_welcome(message):
    msg = message
    bot.reply_to(
        message,
        f"""\
Hey {message.from_user.first_name}!
I've collected your chat id and i'm going to use it to send status updates.\
""",
    )
    global chat_id
    chat_id = message.chat.id
    with open('callid.txt', 'w') as f:
        f.write(str(chat_id))


# Handle all other messages with content_type 'text' (content_types defaults to ['text'])
@bot.message_handler(commands=["status"])
def status_message(message):
    bot.send_message(message.chat.id, message.text)
    print(chat_id)


@bot.message_handler(commands=["track"])
def track(message):
    markup = quick_markup(
        {"model": {"callback_data": "model"}, "comp": {"callback_data": "comp"}}
    )
    reply = bot.send_message(message.chat.id, "Выбрать модель", reply_markup=markup)


@bot.callback_query_handler(func=lambda message: True)
def send_status(call):
    bot.answer_callback_query(call.id, text="Метод выбран")
    if call.data == "model":
        message, image = computing.print_model()
    elif call.data == "comp":
        message, image = computing.print_comp()
    bot.send_message(call.id, message)
    bot.send_photo(call.id, image)


if __name__ == "__main__":
    bot.infinity_polling()
